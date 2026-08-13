import unittest
from unittest.mock import patch

from sqlalchemy.exc import OperationalError as AlchemyOperationalError

from metadome.domain.repositories import GeneRepository
from metadome.domain.error import RecoverableError


class TestRepositories(unittest.TestCase):
    @patch('metadome.domain.repositories.db._make_scoped_session')
    def test_session_always_removed(self, mock_make_scoped_session):
        class FailSession:
            def __init__(self):
                self.removed = False

            def query(self, id_):
                raise Exception("test fail")

            def remove(self):
                self.removed = True

        class Allable:
            def all(self):
                return []

        class Filterable:
            def filter(self, id_):
                return Allable()

        class SuccessSession:
            def __init__(self):
                self.removed = False

            def query(self, id_):
                return Filterable()

            def remove(self):
                self.removed = True

        session = FailSession()
        mock_make_scoped_session.return_value = session

        with self.assertRaises(Exception):
            GeneRepository.retrieve_all_transcript_ids_with_mappings()
        self.assertTrue(session.removed)

        session = SuccessSession()
        mock_make_scoped_session.return_value = session

        GeneRepository.retrieve_all_transcript_ids_with_mappings()
        self.assertTrue(session.removed)

    @patch('metadome.domain.repositories.db._make_scoped_session')
    @patch('metadome.domain.repositories._log.error')
    def test_logs_error(self, mock_log_error, mock_make_scoped_session):
        error_message = "test fail"

        class FailSession:
            def query(self, id_):
                raise Exception(error_message)

            def remove(self):
                pass

        mock_make_scoped_session.return_value = FailSession()

        with self.assertRaises(Exception) as context:
            GeneRepository.retrieve_all_transcript_ids_with_mappings()

        self.assertEqual(str(context.exception), error_message)
        mock_log_error.assert_called()

    @patch('metadome.domain.repositories.db._make_scoped_session')
    def test_raises_recoverable_error(self, mock_make_scoped_session):
        class FailSession:
            def query(self, id_):
                raise AlchemyOperationalError("SELECT 1", {}, Exception("test fail"))

            def remove(self):
                pass

        mock_make_scoped_session.return_value = FailSession()

        with self.assertRaises(RecoverableError):
            GeneRepository.retrieve_all_transcript_ids_with_mappings()


if __name__ == "__main__":
    unittest.main()
