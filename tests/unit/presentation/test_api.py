import json
import logging
import unittest
from unittest.mock import patch

from flask import url_for

from metadome.factory import create_app


_log = logging.getLogger(__name__)


class TestApi(unittest.TestCase):
    @patch("metadome.database.db")
    def setUp(self, mock_db):
        class FakeDB:
            def init_app(self, app):
                pass

            def create_all(self):
                pass

        mock_db.return_value = FakeDB()

        self.server = create_app({
            'SERVER_NAME': "test-server",
            'TESTING': True,
            'SECRET_KEY': 'testing',
        })
        self.client = self.server.test_client()

    @patch("metadome.tasks.retrieve_metadomain_annotation")
    def test_get_metadomain_annotation(self, mock_retrieve):
        mock_retrieve.return_value = {}

        input_ = {
            'transcript_id': 'test',
            'protein_position': 1,
            'requested_domains': {}
        }

        with self.server.app_context():
            response = self.client.post(
                url_for('api.get_metadomain_annotation'),
                data=json.dumps(input_),
                content_type="application/json"
            )

        self.assertEqual(response.status_code, 200)


if __name__ == "__main__":
    unittest.main()
