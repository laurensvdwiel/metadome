import unittest
from unittest.mock import patch

from metadome.factory import create_app, make_celery


class TestTasks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        flask_app = create_app({
            'TESTING': True,
            'CELERYD_CONCURRENCY': 0,
        })
        make_celery(flask_app)

    @patch("metadome.tasks.analyse_transcript")
    def test_update_sent_state(self, mock_analyse):
        mock_analyse.return_value = {'test': "ok"}

        from metadome.tasks import retrieve_prebuild_visualization
        result = retrieve_prebuild_visualization.delay('test_id')

        self.assertEqual(result.status, 'SENT')


if __name__ == "__main__":
    unittest.main()
