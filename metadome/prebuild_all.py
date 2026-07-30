import argparse
import logging
import os
import resource
import time
import traceback
from datetime import datetime, timezone
from multiprocessing import Pool

from metadome.application import app
from metadome.database import db
from metadome.tasks import analyse_transcript
from metadome.controllers.job import (
    store_visualization, store_error,
    get_visualization_path, get_visualization_error_path,
)
from metadome.domain.repositories import GeneRepository

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
logging.getLogger('sqlalchemy.engine').setLevel(logging.WARNING)  # SQLALCHEMY_ECHO=True is very noisy
_log = logging.getLogger(__name__)

_LIST_CAP = 25  # how many example ids to print per category


def _fmt_ts(epoch):
    if epoch is None:
        return "never"
    return datetime.fromtimestamp(epoch, tz=timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')


def _rusage_snapshot():
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime
    peak_rss_kb = max(s.ru_maxrss, c.ru_maxrss)  # Linux reports kilobytes
    return cpu, peak_rss_kb


def _log_run_stats(t0, cpu0, ok, failed):
    wall = time.monotonic() - t0
    cpu1, peak_rss_kb = _rusage_snapshot()
    cpu = cpu1 - cpu0
    done = ok + failed
    _log.info("done: %d ok, %d failed", ok, failed)
    _log.info(
        "stats: wall=%.1fs | %.2f items/s | %.1fs/item avg | "
        "cpu=%.0fs (~%.1f cores busy) | peak_rss(largest proc)=%.0f MB",
        wall,
        (done / wall) if wall > 0 else 0.0,
        (wall / done) if done > 0 else 0.0,
        cpu,
        (cpu / wall) if wall > 0 else 0.0,
        peak_rss_kb / 1024.0,
    )

def _transcript_ids(genome_build):
    return sorted({
        t.gencode_transcription_id
        for t in GeneRepository.retrieve_all_transcript_ids_with_mappings_for_genome_build(genome_build)
    })


def _state(transcript_id, genome_build):
    """Read-only classification: ('built', mtime) | ('error', mtime) | ('missing', None)."""
    vpath = get_visualization_path(transcript_id, genome_build)
    if os.path.isfile(vpath):
        return 'built', os.path.getmtime(vpath)
    epath = get_visualization_error_path(transcript_id, genome_build)
    if os.path.isfile(epath):
        return 'error', os.path.getmtime(epath)
    return 'missing', None


def audit(genome_builds, stale_seconds=None):
    """Stat every expected (transcript, build); return per-build stats. Never writes."""
    now = time.time()
    per_build = {}
    for gb in genome_builds:
        tids = _transcript_ids(gb)
        built = errored = 0
        newest = oldest = None
        missing, error_ids, stale = [], [], []
        for tid in tids:
            state, mtime = _state(tid, gb)
            if state == 'built':
                built += 1
                newest = mtime if newest is None else max(newest, mtime)
                oldest = mtime if oldest is None else min(oldest, mtime)
                if stale_seconds is not None and (now - mtime) > stale_seconds:
                    stale.append(tid)
            elif state == 'error':
                errored += 1
                error_ids.append(tid)
            else:
                missing.append(tid)
        per_build[gb] = {
            'expected': len(tids), 'built': built, 'errored': errored,
            'missing': missing, 'error_ids': error_ids, 'stale': stale,
            'newest': newest, 'oldest': oldest,
        }
    return per_build


def print_report(per_build):
    _log.info("========== PREBUILD AUDIT ==========")
    for gb, s in per_build.items():
        _log.info("--- %s ---", gb)
        _log.info("  expected transcripts : %d", s['expected'])
        _log.info("  built                : %d", s['built'])
        _log.info("  missing              : %d", len(s['missing']))
        _log.info("  errored              : %d", s['errored'])
        _log.info("  oldest build         : %s", _fmt_ts(s['oldest']))
        _log.info("  newest build         : %s", _fmt_ts(s['newest']))
        if s['stale']:
            _log.info("  stale (older cutoff) : %d", len(s['stale']))
        for label, ids in (('missing', s['missing']), ('errored', s['error_ids']), ('stale', s['stale'])):
            if ids:
                shown = ', '.join(ids[:_LIST_CAP])
                more = '' if len(ids) <= _LIST_CAP else f" (+{len(ids) - _LIST_CAP} more)"
                _log.info("  %s ids: %s%s", label, shown, more)
    _log.info("====================================")


def collect_work(genome_builds, overwrite, stale_seconds):
    """Pending = missing + errored (+ stale, + everything if overwrite)."""
    now = time.time()
    work = []
    for gb in genome_builds:
        for tid in _transcript_ids(gb):
            if overwrite:
                work.append((tid, gb))
                continue
            state, mtime = _state(tid, gb)
            if state in ('missing', 'error'):
                work.append((tid, gb))
            elif state == 'built' and stale_seconds is not None and (now - mtime) > stale_seconds:
                work.append((tid, gb))
    return work


def _prebuild_one(item):
    transcript_id, genome_build = item
    try:
        result = analyse_transcript(transcript_id, genome_build)
        if 'error' in result:
            store_error(transcript_id, genome_build, result['error'])
            return (genome_build, transcript_id, 'ERROR', result['error'])
        store_visualization(transcript_id, genome_build, result)  # LockFile-protected
        return (genome_build, transcript_id, 'OK', None)
    except Exception:
        tb = traceback.format_exc()
        store_error(transcript_id, genome_build, tb)
        return (genome_build, transcript_id, 'FAIL', tb.splitlines()[-1])


def _init_worker():
    # Forked workers must not reuse the parent's DB connections.
    app.app_context().push()
    db.engine.dispose()


def _run(work, workers):
    total, ok, failed = len(work), 0, 0
    _log.info("prebuilding %d (transcript, genome_build) pairs with %d worker(s)", total, workers)

    def _tally(i, res):
        nonlocal ok, failed
        gb, tid, status, msg = res
        if status == 'OK':
            ok += 1
        else:
            failed += 1
            _log.error("%s/%s: %s %s", gb, tid, status, msg)
        if i % 100 == 0:
            _log.info("progress: %d/%d (ok=%d failed=%d)", i, total, ok, failed)

    t0 = time.monotonic()
    cpu0, _ = _rusage_snapshot()
    try:
        if workers <= 1:
            with app.app_context():
                for i, item in enumerate(work, 1):
                    _tally(i, _prebuild_one(item))
        else:
            with Pool(processes=workers, initializer=_init_worker) as pool:
                try:
                    for i, res in enumerate(pool.imap_unordered(_prebuild_one, work), 1):
                        _tally(i, res)
                except KeyboardInterrupt:
                    _log.warning("interrupt received — terminating workers...")
                    pool.terminate()
                    pool.join()
                    raise
    except KeyboardInterrupt:
        _log.warning("STOPPED EARLY (already-built files are kept).")
    finally:
        _log_run_stats(t0, cpu0, ok, failed)


def main():
    p = argparse.ArgumentParser(description="Audit and prebuild transcript visualizations for both genome builds.")
    p.add_argument('--report', action='store_true', help='Audit only (read-only); do not build anything.')
    p.add_argument('--workers', type=int, default=1, help='Parallel worker processes (default 1).')
    p.add_argument('--limit', type=int, default=None, help='Build at most N pending items, then stop.')
    p.add_argument('--genome-build', action='append', default=None,
                   help='Restrict to a build (repeatable). Default: all builds in the DB.')
    p.add_argument('--overwrite', action='store_true', help='Rebuild every transcript, even if already built.')
    p.add_argument('--stale-days', type=float, default=None,
                   help='Treat visualizations older than N days as pending (rebuild them).')
    args = p.parse_args()

    stale_seconds = args.stale_days * 86400 if args.stale_days is not None else None

    with app.app_context():
        genome_builds = args.genome_build or GeneRepository.retrieve_all_genome_builds_from_db()
        _log.info("genome builds: %s", genome_builds)

        # Always audit first — answers "when were prebuilds last made?" and checks every transcript/build.
        print_report(audit(genome_builds, stale_seconds))

        if args.report:
            return

        work = collect_work(genome_builds, args.overwrite, stale_seconds)
        if args.limit is not None:
            _log.info("limiting to first %d of %d pending", args.limit, len(work))
            work = work[:args.limit]

    if not work:
        _log.info("nothing to build — all transcripts up to date.")
        return
    _run(work, args.workers)


if __name__ == '__main__':
    main()