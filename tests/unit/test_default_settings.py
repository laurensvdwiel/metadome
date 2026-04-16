import unittest

import metadome.default_settings as settings


COVERED_SETTINGS = {
    "SECRET_KEY",
    "GOOGLE_SITE_VERIFICATION",
    "GOOGLE_ANALYTICS_ID",
    "RECAPTCHA_SITE_KEY",
    "RECAPTCHA_SECRET_KEY",
    "RECAPTCHA_THRESHOLD",
    "DEBUG",
    "SQLALCHEMY_RECORD_QUERIES",
    "SQLALCHEMY_TRACK_MODIFICATIONS",
    "SQLALCHEMY_DATABASE_URI",
    "SQLALCHEMY_ECHO",
    "SQLALCHEMY_POOL_TIMEOUT",
    "SQLALCHEMY_ENGINE_OPTIONS",
    "CELERY_BROKER_URL",
    "CELERY_RESULT_BACKEND",
    "CELERY_TRACK_STARTED",
    "CELERY_TASK_SERIALIZER",
    "CELERY_RESULT_SERIALIZER",
    "CELERY_ACCEPT_CONTENT",
    "CACHE_TYPE",
    "CACHE_REDIS_URL",
    "CACHE_DEFAULT_TIMEOUT",
    "REPOSITORY_CACHE_TIMEOUT",
    "ALLELE_FREQUENCY_CUTOFF",
    "SLIDING_WINDOW_SIZE",
    "DEBUG_TB_ENABLED",
    "MAIL_SERVER",
    "MAIL_PORT",
    "MAIL_USE_TLS",
    "MAIL_USE_SSL",
    "SUPPORT_EMAIL",
    "ISSUES_EMAIL",
    "MAIL_DEFAULT_SENDER",
    "ERROR_EMAIL_NOTIFICATION_WINDOW",
    "DATA_DIR",
    "GENE_NAMES_FILE",
    "GENOME_BUILDS_FILE",
    "BLASTP_EXECUTABLE",
    "CLUSTALW_EXECUTABLE",
    "HMMFETCH_EXECUTABLE",
    "HMMLOGO_EXECUTABLE",
    "HMMALIGN_EXECUTABLE",
    "HMMEMIT_EXECUTABLE",
    "HMMSTAT_EXECUTABLE",
    "GENCODE_HG_ANNOTATION_FILE_GTF",
    "GENCODE_HG_ANNOTATION_FILE_GFF3",
    "GENCODE_HG_TRANSCRIPTION_FILE",
    "GENCODE_HG_TRANSLATION_FILE",
    "GENCODE_REFSEQ_FILE",
    "GENCODE_SWISSPROT_FILE",
    "GENCODE_BASIC_FILE",
    "INTERPROSCAN_DOCKER_IMAGE",
    "INTERPROSCAN_DOCKER_VOLUME",
    "INTERPROSCAN_EXECUTABLE",
    "INTERPROSCAN_TEMP_DIR",
    "INTERPROSCAN_DOMAIN_DATABASES",
    "UNIPROT_MAX_BLAST_RESULTS",
    "UNIPROT_DIR",
    "UNIPROT_SPROT_CANONICAL_AND_ISOFORM",
    "UNIPROT_SPROT_ISOFORM",
    "UNIPROT_SPROT_CANONICAL",
    "UNIPROT_SPROT_SPECIES_FILTER",
    "METADOMAIN_DIR",
    "RECONSTRUCT_METADOMAINS",
    "METADOMAIN_ALIGNMENT_FILE_NAME",
    "METADOMAIN_MAPPING_FILE_NAME",
    "METADOMAIN_DETAILS_FILE_NAME",
    "METADOMAIN_SNV_ANNOTATION_FILE_NAME",
    "PRE_BUILD_VISUALIZATION_DIR",
    "PRE_BUILD_VISUALIZATION_FILE_NAME",
    "PRE_BUILD_VISUALIZATION_TASK_FILE_NAME",
    "PRE_BUILD_VISUALIZATION_ERROR_FILE_NAME",
    "PFAM_DIR",
    "PFAM_ALIGNMENT_DIR",
    "PFAM_HMM_DAT",
    "PFAM_HMM",
    "GNOMAD_DIR",
    "GNOMAD_GRCH37_VCF_FILE",
    "GNOMAD_GRCH38_VCF_FILE",
    "GNOMAD_ACCEPTED_FILTERS",
    "CLINVAR_DIR",
    "CLINVAR_GRCH37_VCF_FILE",
    "CLINVAR_GRCH38_VCF_FILE",
    "CLINVAR_CONSIDERED_CLINSIG",
}


def _public_setting_names():
    return {
        name for name in dir(settings)
        if name.isupper()
    }


def _assert_non_empty_string(testcase, value, name):
    testcase.assertIsInstance(value, str, f"{name} must be a string")
    testcase.assertTrue(value.strip(), f"{name} must not be empty")


def _assert_optional_string(testcase, value, name):
    testcase.assertTrue(
        value is None or (isinstance(value, str) and value.strip()),
        f"{name} must be None or a non-empty string",
    )


def _assert_positive_int(testcase, value, name):
    testcase.assertIsInstance(value, int, f"{name} must be an int")
    testcase.assertGreater(value, 0, f"{name} must be > 0")


def _assert_non_negative_number(testcase, value, name):
    testcase.assertIsInstance(value, (int, float), f"{name} must be numeric")
    testcase.assertGreaterEqual(value, 0, f"{name} must be >= 0")


class TestDefaultSettings(unittest.TestCase):
    def test_all_public_settings_are_explicitly_covered(self):
        current_settings = _public_setting_names()

        missing_from_tests = sorted(current_settings - COVERED_SETTINGS)
        stale_in_tests = sorted(COVERED_SETTINGS - current_settings)

        self.assertEqual(
            [],
            missing_from_tests,
            msg=(
                "New settings were added to metadome.default_settings but are not yet "
                "covered by tests: "
                + ", ".join(missing_from_tests)
            ),
        )
        self.assertEqual(
            [],
            stale_in_tests,
            msg=(
                "Tests still reference settings that no longer exist in "
                "metadome.default_settings: "
                + ", ".join(stale_in_tests)
            ),
        )

    def test_core_boolean_settings(self):
        self.assertIs(settings.DEBUG, False)
        self.assertIs(settings.SQLALCHEMY_TRACK_MODIFICATIONS, False)
        self.assertIs(settings.SQLALCHEMY_RECORD_QUERIES, settings.DEBUG)
        self.assertIs(settings.DEBUG_TB_ENABLED, settings.DEBUG)
        self.assertIsInstance(settings.SQLALCHEMY_ECHO, bool)
        self.assertIsInstance(settings.CELERY_TRACK_STARTED, bool)
        self.assertIsInstance(settings.MAIL_USE_TLS, bool)
        self.assertIsInstance(settings.MAIL_USE_SSL, bool)
        self.assertIsInstance(settings.RECONSTRUCT_METADOMAINS, bool)

    def test_core_numeric_settings(self):
        self.assertIsInstance(settings.RECAPTCHA_THRESHOLD, (int, float))
        self.assertGreaterEqual(settings.RECAPTCHA_THRESHOLD, 0.0)
        self.assertLessEqual(settings.RECAPTCHA_THRESHOLD, 1.0)

        _assert_positive_int(self, settings.SQLALCHEMY_POOL_TIMEOUT, "SQLALCHEMY_POOL_TIMEOUT")
        _assert_positive_int(self, settings.CACHE_DEFAULT_TIMEOUT, "CACHE_DEFAULT_TIMEOUT")
        _assert_positive_int(self, settings.REPOSITORY_CACHE_TIMEOUT, "REPOSITORY_CACHE_TIMEOUT")
        _assert_positive_int(self, settings.SLIDING_WINDOW_SIZE, "SLIDING_WINDOW_SIZE")
        _assert_positive_int(self, settings.ERROR_EMAIL_NOTIFICATION_WINDOW, "ERROR_EMAIL_NOTIFICATION_WINDOW")
        _assert_positive_int(self, settings.MAIL_PORT, "MAIL_PORT")
        _assert_positive_int(self, settings.UNIPROT_MAX_BLAST_RESULTS, "UNIPROT_MAX_BLAST_RESULTS")

        _assert_non_negative_number(self, settings.ALLELE_FREQUENCY_CUTOFF, "ALLELE_FREQUENCY_CUTOFF")

    def test_core_string_and_optional_string_settings(self):
        _assert_non_empty_string(self, settings.SQLALCHEMY_DATABASE_URI, "SQLALCHEMY_DATABASE_URI")
        self.assertTrue(
            settings.SQLALCHEMY_DATABASE_URI.startswith("postgresql+psycopg://"),
            "SQLALCHEMY_DATABASE_URI should use the PostgreSQL scheme",
        )

        _assert_non_empty_string(self, settings.CELERY_BROKER_URL, "CELERY_BROKER_URL")
        _assert_non_empty_string(self, settings.CELERY_RESULT_BACKEND, "CELERY_RESULT_BACKEND")
        _assert_non_empty_string(self, settings.CACHE_TYPE, "CACHE_TYPE")
        _assert_non_empty_string(self, settings.CACHE_REDIS_URL, "CACHE_REDIS_URL")
        _assert_non_empty_string(self, settings.DATA_DIR, "DATA_DIR")

        _assert_optional_string(self, settings.MAIL_SERVER, "MAIL_SERVER")
        _assert_optional_string(self, settings.SUPPORT_EMAIL, "SUPPORT_EMAIL")
        _assert_optional_string(self, settings.ISSUES_EMAIL, "ISSUES_EMAIL")
        _assert_optional_string(self, settings.MAIL_DEFAULT_SENDER, "MAIL_DEFAULT_SENDER")
        _assert_optional_string(self, settings.GOOGLE_SITE_VERIFICATION, "GOOGLE_SITE_VERIFICATION")
        _assert_optional_string(self, settings.GOOGLE_ANALYTICS_ID, "GOOGLE_ANALYTICS_ID")
        _assert_optional_string(self, settings.RECAPTCHA_SITE_KEY, "RECAPTCHA_SITE_KEY")
        _assert_optional_string(self, settings.RECAPTCHA_SECRET_KEY, "RECAPTCHA_SECRET_KEY")

        self.assertIsInstance(settings.SECRET_KEY, str)
        self.assertTrue(len(settings.SECRET_KEY) > 0, "SECRET_KEY must not be empty")

    def test_sqlalchemy_engine_options(self):
        engine_options = settings.SQLALCHEMY_ENGINE_OPTIONS

        self.assertIsInstance(engine_options, dict)
        self.assertEqual(
            {"pool_pre_ping", "pool_recycle", "pool_timeout", "connect_args"},
            set(engine_options.keys()),
        )
        self.assertIs(engine_options["pool_pre_ping"], True)
        _assert_positive_int(self, engine_options["pool_recycle"], "SQLALCHEMY_ENGINE_OPTIONS['pool_recycle']")
        self.assertEqual(engine_options["pool_timeout"], settings.SQLALCHEMY_POOL_TIMEOUT)

        connect_args = engine_options["connect_args"]
        self.assertIsInstance(connect_args, dict)
        self.assertEqual(
            {
                "application_name",
                "connect_timeout",
                "keepalives",
                "keepalives_idle",
                "keepalives_interval",
                "keepalives_count",
            },
            set(connect_args.keys()),
        )
        self.assertEqual(connect_args["application_name"], "metadome")
        _assert_positive_int(self, connect_args["connect_timeout"], "connect_timeout")
        self.assertEqual(connect_args["keepalives"], 1)
        _assert_positive_int(self, connect_args["keepalives_idle"], "keepalives_idle")
        _assert_positive_int(self, connect_args["keepalives_interval"], "keepalives_interval")
        _assert_positive_int(self, connect_args["keepalives_count"], "keepalives_count")

    def test_celery_settings(self):
        self.assertEqual(settings.CELERY_TASK_SERIALIZER, "pickle")
        self.assertEqual(settings.CELERY_RESULT_SERIALIZER, "pickle")
        self.assertIsInstance(settings.CELERY_ACCEPT_CONTENT, list)
        self.assertIn("pickle", settings.CELERY_ACCEPT_CONTENT)

    def test_cache_settings(self):
        self.assertEqual(settings.CACHE_TYPE, "RedisCache")
        self.assertTrue(
            settings.CACHE_REDIS_URL.startswith("redis://"),
            "CACHE_REDIS_URL should use the Redis URL scheme",
        )
        self.assertGreaterEqual(settings.REPOSITORY_CACHE_TIMEOUT, settings.CACHE_DEFAULT_TIMEOUT)

    def test_directory_path_conventions(self):
        directory_settings = {
            "DATA_DIR": settings.DATA_DIR,
            "UNIPROT_DIR": settings.UNIPROT_DIR,
            "METADOMAIN_DIR": settings.METADOMAIN_DIR,
            "PRE_BUILD_VISUALIZATION_DIR": settings.PRE_BUILD_VISUALIZATION_DIR,
            "PFAM_DIR": settings.PFAM_DIR,
            "PFAM_ALIGNMENT_DIR": settings.PFAM_ALIGNMENT_DIR,
            "GNOMAD_DIR": settings.GNOMAD_DIR,
            "CLINVAR_DIR": settings.CLINVAR_DIR,
            "INTERPROSCAN_TEMP_DIR": settings.INTERPROSCAN_TEMP_DIR,
        }

        for name, value in directory_settings.items():
            _assert_non_empty_string(self, value, name)

        self.assertTrue(settings.DATA_DIR.endswith("/"))
        self.assertTrue(settings.UNIPROT_DIR.startswith(settings.DATA_DIR))
        self.assertTrue(settings.METADOMAIN_DIR.startswith(settings.DATA_DIR))
        self.assertTrue(settings.PRE_BUILD_VISUALIZATION_DIR.startswith(settings.DATA_DIR))
        self.assertTrue(settings.PFAM_DIR.startswith(settings.DATA_DIR))
        self.assertTrue(settings.PFAM_ALIGNMENT_DIR.startswith(settings.PFAM_DIR))
        self.assertTrue(settings.GNOMAD_DIR.startswith(settings.DATA_DIR))
        self.assertTrue(settings.CLINVAR_DIR.startswith(settings.DATA_DIR))
        self.assertTrue(settings.INTERPROSCAN_TEMP_DIR.startswith("/"))

    def test_file_path_conventions(self):
        file_settings = {
            "GENE_NAMES_FILE": settings.GENE_NAMES_FILE,
            "GENOME_BUILDS_FILE": settings.GENOME_BUILDS_FILE,
            "BLASTP_EXECUTABLE": settings.BLASTP_EXECUTABLE,
            "CLUSTALW_EXECUTABLE": settings.CLUSTALW_EXECUTABLE,
            "HMMFETCH_EXECUTABLE": settings.HMMFETCH_EXECUTABLE,
            "HMMLOGO_EXECUTABLE": settings.HMMLOGO_EXECUTABLE,
            "HMMALIGN_EXECUTABLE": settings.HMMALIGN_EXECUTABLE,
            "HMMEMIT_EXECUTABLE": settings.HMMEMIT_EXECUTABLE,
            "HMMSTAT_EXECUTABLE": settings.HMMSTAT_EXECUTABLE,
            "GENCODE_HG_ANNOTATION_FILE_GTF": settings.GENCODE_HG_ANNOTATION_FILE_GTF,
            "GENCODE_HG_ANNOTATION_FILE_GFF3": settings.GENCODE_HG_ANNOTATION_FILE_GFF3,
            "GENCODE_HG_TRANSCRIPTION_FILE": settings.GENCODE_HG_TRANSCRIPTION_FILE,
            "GENCODE_HG_TRANSLATION_FILE": settings.GENCODE_HG_TRANSLATION_FILE,
            "GENCODE_REFSEQ_FILE": settings.GENCODE_REFSEQ_FILE,
            "GENCODE_SWISSPROT_FILE": settings.GENCODE_SWISSPROT_FILE,
            "GENCODE_BASIC_FILE": settings.GENCODE_BASIC_FILE,
            "UNIPROT_SPROT_CANONICAL_AND_ISOFORM": settings.UNIPROT_SPROT_CANONICAL_AND_ISOFORM,
            "UNIPROT_SPROT_ISOFORM": settings.UNIPROT_SPROT_ISOFORM,
            "UNIPROT_SPROT_CANONICAL": settings.UNIPROT_SPROT_CANONICAL,
            "PFAM_HMM_DAT": settings.PFAM_HMM_DAT,
            "PFAM_HMM": settings.PFAM_HMM,
            "GNOMAD_GRCH37_VCF_FILE": settings.GNOMAD_GRCH37_VCF_FILE,
            "GNOMAD_GRCH38_VCF_FILE": settings.GNOMAD_GRCH38_VCF_FILE,
            "CLINVAR_GRCH37_VCF_FILE": settings.CLINVAR_GRCH37_VCF_FILE,
            "CLINVAR_GRCH38_VCF_FILE": settings.CLINVAR_GRCH38_VCF_FILE,
        }

        for name, value in file_settings.items():
            _assert_non_empty_string(self, value, name)
            self.assertIn("/", value, f"{name} should look like a path")

        self.assertTrue(settings.GENE_NAMES_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENOME_BUILDS_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_HG_ANNOTATION_FILE_GTF.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_HG_ANNOTATION_FILE_GFF3.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_HG_TRANSCRIPTION_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_HG_TRANSLATION_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_REFSEQ_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_SWISSPROT_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.GENCODE_BASIC_FILE.startswith(settings.DATA_DIR))
        self.assertTrue(settings.UNIPROT_SPROT_CANONICAL_AND_ISOFORM.startswith(settings.UNIPROT_DIR))
        self.assertTrue(settings.UNIPROT_SPROT_ISOFORM.startswith(settings.UNIPROT_DIR))
        self.assertTrue(settings.UNIPROT_SPROT_CANONICAL.startswith(settings.UNIPROT_DIR))
        self.assertTrue(settings.PFAM_HMM_DAT.startswith(settings.PFAM_DIR))
        self.assertTrue(settings.PFAM_HMM.startswith(settings.PFAM_DIR))
        self.assertTrue(settings.GNOMAD_GRCH37_VCF_FILE.startswith(settings.GNOMAD_DIR))
        self.assertTrue(settings.GNOMAD_GRCH38_VCF_FILE.startswith(settings.GNOMAD_DIR))
        self.assertTrue(settings.CLINVAR_GRCH37_VCF_FILE.startswith(settings.CLINVAR_DIR))
        self.assertTrue(settings.CLINVAR_GRCH38_VCF_FILE.startswith(settings.CLINVAR_DIR))

    def test_file_name_settings(self):
        file_name_settings = {
            "METADOMAIN_ALIGNMENT_FILE_NAME": settings.METADOMAIN_ALIGNMENT_FILE_NAME,
            "METADOMAIN_MAPPING_FILE_NAME": settings.METADOMAIN_MAPPING_FILE_NAME,
            "METADOMAIN_DETAILS_FILE_NAME": settings.METADOMAIN_DETAILS_FILE_NAME,
            "METADOMAIN_SNV_ANNOTATION_FILE_NAME": settings.METADOMAIN_SNV_ANNOTATION_FILE_NAME,
            "PRE_BUILD_VISUALIZATION_FILE_NAME": settings.PRE_BUILD_VISUALIZATION_FILE_NAME,
            "PRE_BUILD_VISUALIZATION_TASK_FILE_NAME": settings.PRE_BUILD_VISUALIZATION_TASK_FILE_NAME,
            "PRE_BUILD_VISUALIZATION_ERROR_FILE_NAME": settings.PRE_BUILD_VISUALIZATION_ERROR_FILE_NAME,
        }

        for name, value in file_name_settings.items():
            _assert_non_empty_string(self, value, name)
            self.assertNotIn("/", value, f"{name} should be a file name, not a full path")

        self.assertTrue(settings.METADOMAIN_DETAILS_FILE_NAME.endswith(".json"))
        self.assertTrue(settings.PRE_BUILD_VISUALIZATION_FILE_NAME.endswith(".json"))

    def test_domain_specific_collection_settings(self):
        self.assertEqual(settings.INTERPROSCAN_DOMAIN_DATABASES, "Pfam")
        self.assertEqual(settings.UNIPROT_SPROT_SPECIES_FILTER, "HUMAN")

        self.assertIsInstance(settings.GNOMAD_ACCEPTED_FILTERS, list)
        self.assertGreater(len(settings.GNOMAD_ACCEPTED_FILTERS), 0)
        self.assertTrue(all(isinstance(x, str) and x for x in settings.GNOMAD_ACCEPTED_FILTERS))

        self.assertIsInstance(settings.CLINVAR_CONSIDERED_CLINSIG, list)
        self.assertGreater(len(settings.CLINVAR_CONSIDERED_CLINSIG), 0)
        self.assertTrue(all(isinstance(x, str) and x for x in settings.CLINVAR_CONSIDERED_CLINSIG))

    def test_interproscan_settings(self):
        _assert_non_empty_string(self, settings.INTERPROSCAN_DOCKER_IMAGE, "INTERPROSCAN_DOCKER_IMAGE")
        _assert_non_empty_string(self, settings.INTERPROSCAN_DOCKER_VOLUME, "INTERPROSCAN_DOCKER_VOLUME")
        _assert_non_empty_string(self, settings.INTERPROSCAN_EXECUTABLE, "INTERPROSCAN_EXECUTABLE")

        self.assertIn(":", settings.INTERPROSCAN_DOCKER_IMAGE)
        self.assertTrue(settings.INTERPROSCAN_EXECUTABLE.endswith(".sh"))


if __name__ == "__main__":
    unittest.main()