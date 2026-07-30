import os
import logging
import json

import tempfile
from lockfile import LockFile

from celery.result import AsyncResult
from flask import current_app as flask_app


_log = logging.getLogger(__name__)


def get_visualization_path(transcript_id, genome_build):
    return os.path.join(_get_visualization_dir_path(transcript_id, genome_build),
                        flask_app.config['PRE_BUILD_VISUALIZATION_FILE_NAME'])


def _get_visualization_dir_path(transcript_id, genome_build):
    return os.path.join(flask_app.config['PRE_BUILD_VISUALIZATION_DIR'], genome_build, transcript_id)


def _get_visualization_task_path(transcript_id, genome_build):
    return os.path.join(_get_visualization_dir_path(transcript_id, genome_build),
                        flask_app.config['PRE_BUILD_VISUALIZATION_TASK_FILE_NAME'])


def get_visualization_error_path(transcript_id, genome_build):
    return os.path.join(_get_visualization_dir_path(transcript_id, genome_build),
                        flask_app.config['PRE_BUILD_VISUALIZATION_ERROR_FILE_NAME'])


def _get_lock_for(transcript_id, genome_build):
    lock_dir_path = _get_visualization_dir_path(transcript_id, genome_build)

    if not os.path.isdir(lock_dir_path):
        os.makedirs(lock_dir_path)

    return LockFile(lock_dir_path)


def _cleanup_visualisation_if_needed(transcript_id, genome_build):
    visualization_path = get_visualization_path(transcript_id, genome_build)
    task_path = _get_visualization_task_path(transcript_id, genome_build)

    if os.path.isfile(visualization_path) and os.path.isfile(task_path):
        # Not needed anymore
        os.remove(task_path)


def create_visualization_job_if_needed(transcript_id, genome_build):
    visualization_path = get_visualization_path(transcript_id, genome_build)
    task_path = _get_visualization_task_path(transcript_id, genome_build)
    error_path = get_visualization_error_path(transcript_id, genome_build)

    with _get_lock_for(transcript_id, genome_build):
        if os.path.isfile(error_path):
            # It has failed before, try this job again.
            os.remove(error_path)
            if os.path.isfile(task_path):
                os.remove(task_path)

        elif os.path.isfile(task_path):
            with open(task_path, 'r') as f:
                task_id = f.read()

            result = AsyncResult(task_id)
            if result.status == 'PENDING':  # PENDING means it's just not in the backend
                os.remove(task_path)
            else:
                _log.info("visualization job for genome_build {} transcript {} is already submitted as task {}"
                          .format(genome_build, transcript_id, task_id))
                return

        if os.path.isfile(visualization_path):
            _log.info("visualization file for genome_build {} transcript {} already exists"
                      .format(genome_build, transcript_id))
        else:
            from metadome.tasks import create_prebuild_visualization
            result = create_prebuild_visualization.delay(transcript_id, genome_build)

            with open(task_path, 'w') as f:
                f.write(result.task_id)

            # From here on, the task itself will handle the creation of result and error files.


def get_visualization_status(transcript_id, genome_build):
    visualization_path = get_visualization_path(transcript_id, genome_build)
    error_path = get_visualization_error_path(transcript_id, genome_build)
    task_path = _get_visualization_task_path(transcript_id, genome_build)

    with _get_lock_for(transcript_id, genome_build):
        _cleanup_visualisation_if_needed(transcript_id, genome_build)

        if os.path.isfile(visualization_path):
            return 'SUCCESS'
        elif os.path.isfile(error_path):
            return 'FAILURE'
        elif os.path.isfile(task_path):
            with open(task_path, 'r') as f:
                task_id = f.read()
                result = AsyncResult(task_id)

                return result.status
        else:
            return 'PENDING'


def store_error(transcript_id, genome_build, traceback):
    error_path = get_visualization_error_path(transcript_id, genome_build)

    with _get_lock_for(transcript_id, genome_build):
        with open(error_path, 'w') as f:
            f.write(traceback)


def retrieve_error(transcript_id, genome_build):
    error_path = get_visualization_error_path(transcript_id, genome_build)

    with _get_lock_for(transcript_id, genome_build):
        if os.path.isfile(error_path):
            with open(error_path, 'r') as f:
                # Get file modification time
                error_timestamp = os.path.getmtime(error_path)

                return f.read(), error_timestamp
        else:
            return 'unknown', None


def store_visualization(transcript_id, genome_build, result):
    visualization_path = get_visualization_path(transcript_id, genome_build)

    with _get_lock_for(transcript_id, genome_build):
        dir_name = os.path.dirname(visualization_path)
        with tempfile.NamedTemporaryFile('w', dir=dir_name, delete=False) as f:
            json.dump(result, f)
            f.flush()
            os.fsync(f.fileno())
            tmp_path = f.name
        os.replace(tmp_path, visualization_path)  # atomic swap


def retrieve_visualization(transcript_id, genome_build):
    visualization_path = get_visualization_path(transcript_id, genome_build)

    with _get_lock_for(transcript_id, genome_build):
        _cleanup_visualisation_if_needed(transcript_id, genome_build)

        if not os.path.isfile(visualization_path):
            raise FileNotFoundError("missing file: {}".format(visualization_path))

        with open(visualization_path, 'r') as f:
            visualization_content = json.load(f)

            return visualization_content
