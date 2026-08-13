from shutil import rmtree
from threading import Thread
from tempfile import mkdtemp
from time import sleep
import logging
import os
import traceback
import unittest
from unittest.mock import patch

from metadome.factory import create_app
from metadome.controllers.job import (
    create_visualization_job_if_needed,
    get_visualization_status,
    retrieve_visualization,
    get_visualization_path,
    store_visualization,
    store_error,
)


_log = logging.getLogger(__name__)


class TestJob(unittest.TestCase):
    @patch("metadome.database.db.create_all")
    def setUp(self, mock_create_all):
        self.temp_dir = mkdtemp()
        mock_create_all.return_value = None
        self.app = create_app()

    def tearDown(self):
        if os.path.isdir(self.temp_dir):
            rmtree(self.temp_dir)

    @patch("metadome.tasks.create_prebuild_visualization.delay")
    @patch("celery.result.AsyncResult")
    @patch("metadome.controllers.job._get_visualization_dir_path")
    def test_run(self, mock_dir_path, mock_result, mock_delay):
        mock_dir_path.return_value = self.temp_dir
        app = self.app

        class ThreadResult(Thread):
            def __init__(self, transcript_id):
                super().__init__()
                self.name = "result_test"
                self.transcript_id = transcript_id
                self.status = "PENDING"

            def run(self):
                self.status = "STARTED"
                try:
                    sleep(5.0)
                    with app.app_context():
                        store_visualization(self.transcript_id, {'id': 'test'})
                    self.status = "SUCCESS"
                except Exception:
                    with app.app_context():
                        store_error(self.transcript_id, traceback.format_exc())
                    self.status = "FAILURE"

        transcript_id = "test_transcript"
        result_thread = ThreadResult(transcript_id)

        mock_result.return_value = result_thread
        mock_delay.return_value = result_thread

        with self.app.app_context():
            create_visualization_job_if_needed(transcript_id)
            result_thread.start()

            while True:
                status = get_visualization_status(transcript_id)
                self.assertNotEqual(status, 'FAILURE')

                if status == 'SUCCESS':
                    break

                sleep(1.0)

            result = retrieve_visualization(transcript_id)

            self.assertGreater(len(result), 0)
            self.assertTrue(os.path.isfile(get_visualization_path(transcript_id)))


if __name__ == "__main__":
    unittest.main()
