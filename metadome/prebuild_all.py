import traceback
import logging
import os
from time import sleep

import metadome.default_settings as settings
from metadome.application import app, celery
from metadome.tasks import create_prebuild_visualization, initialize_metadomain
from metadome.domain.repositories import GeneRepository

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
_log = logging.getLogger(__name__)


def not_created(transcript_id):
    path = os.path.join(settings.PRE_BUILD_VISUALIZATION_DIR,
                        transcript_id,
                        settings.PRE_BUILD_VISUALIZATION_FILE_NAME)
    return not os.path.isfile(path)


def monitor(input_results):
    count_success = 0
    count_failure = 0
    while len(input_results) > 0:
        waiting = {}
        count_running = 0
        for id_ in input_results:

            if input_results[id_].status in ['FAILURE', 'REVOKED']:
                count_failure += 1

                _log.error("{}: {}".format(id_, input_results[id_].traceback))

            elif input_results[id_].status == 'SUCCESS':
                count_success += 1

                try:
                    _log.info("{}: {}".format(id_, input_results[id_].get()))
                except:
                    _log.error(traceback.format_exc())
            else:
                if input_results[id_].status == 'STARTED':
                    count_running += 1

                waiting[id_] = input_results[id_]

        _log.info("failure: {}, success: {}, waiting: {} (of which {} running)"
                  .format(count_failure, count_success,
                          len(waiting), count_running))
        input_results = waiting
        sleep(10)

    _log.debug("all done, stopping..")


_log.debug("submitting metadomain jobs")
results = {domain_id: initialize_metadomain.delay(domain_id)
           for domain_id in os.listdir(settings.METADOMAIN_DIR)}

_log.debug("waiting for results")
try:
    monitor(results)
except:
    _log.debug("revoking all jobs")
    for result in results.values():
        result.revoke()
    raise


_log.debug("getting transcript ids")
with app.app_context():
    transcripts = GeneRepository.retrieve_all_transcript_ids_with_mappings()

    transcript_ids = filter(not_created,
                            [transcript.gencode_transcription_id for transcript in transcripts])

_log.debug("submitting visualization jobs")
results = {transcript_id: create_prebuild_visualization.delay(transcript_id)
           for transcript_id in transcript_ids}

_log.debug("waiting for results")
try:
    monitor(results)
except:
    _log.debug("revoking all jobs")
    for result in results.values():
        result.revoke()
    raise
