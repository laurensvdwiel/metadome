# Credentials
import metadome.flask_app_credentials as credentials
SECRET_KEY = credentials.SECRET_KEY_CRED
GOOGLE_SITE_VERIFICATION = credentials.GOOGLE_SITE_VERIFICATION
GOOGLE_ANALYTICS_ID = credentials.GOOGLE_ANALYTICS_ID
RECAPTCHA_SITE_KEY =  credentials.RECAPTCHA_SITE_KEY
RECAPTCHA_SECRET_KEY = credentials.RECAPTCHA_SECRET_KEY
RECAPTCHA_THRESHOLD = credentials.RECAPTCHA_THRESHOLD

# Flask settings
DEBUG = False

# FLask-SQLAchemy settings
import os
SQLALCHEMY_RECORD_QUERIES = DEBUG # should be false when not debug
SQLALCHEMY_TRACK_MODIFICATIONS = False
SQLALCHEMY_DATABASE_URI = "postgresql+psycopg://"+os.environ["POSTGRES_USER"]+":"+os.environ["POSTGRES_PASSWORD"]+"@"+os.environ["POSTGRES_HOST"]+"/"+os.environ["POSTGRES_DB"]
SQLALCHEMY_ECHO = True
SQLALCHEMY_POOL_TIMEOUT = 10
SQLALCHEMY_ENGINE_OPTIONS = {
    "pool_pre_ping": True,
    "pool_recycle": 1800,
    "pool_timeout": SQLALCHEMY_POOL_TIMEOUT,
    "connect_args": {
        "application_name": "metadome",
        "connect_timeout": 10,
        "keepalives": 1,
        "keepalives_idle": 30,
        "keepalives_interval": 10,
        "keepalives_count": 5,
    },
}

# Flask-Celery settings
CELERY_BROKER_URL='amqp://guest@metadome-dev-rabbitmq-1'
CELERY_RESULT_BACKEND='redis://metadome-dev-redis-1/0'
CELERY_TRACK_STARTED = True
CELERY_TASK_SERIALIZER = 'pickle'
CELERY_RESULT_SERIALIZER='pickle'
CELERY_ACCEPT_CONTENT = ['pickle']
CELERY_TASK_MAX_RETRIES = 50

# Flask-Caching Settings
CACHE_TYPE = 'RedisCache'
CACHE_REDIS_URL = os.environ.get('CACHE_REDIS_URL', 'redis://metadome-dev-redis-1/1')
CACHE_DEFAULT_TIMEOUT = int(os.environ.get('CACHE_DEFAULT_TIMEOUT', '3600'))
REPOSITORY_CACHE_TIMEOUT = int(os.environ.get('REPOSITORY_CACHE_TIMEOUT', '86400'))

# Visualiation specific settings
ALLELE_FREQUENCY_CUTOFF = 0.0
SLIDING_WINDOW_SIZE = 10

# Debug toolbar
DEBUG_TB_ENABLED = DEBUG

# E-mail configuration
MAIL_SERVER = os.environ.get('MAIL_SERVER')
MAIL_PORT = int(os.environ.get('MAIL_PORT', 25))
MAIL_USE_TLS = False
MAIL_USE_SSL = False
# Email addresses from environment
SUPPORT_EMAIL = os.environ.get('SUPPORT_EMAIL')
ISSUES_EMAIL = os.environ.get('ISSUES_EMAIL')
MAIL_DEFAULT_SENDER = os.environ.get('DEFAULT_FROM_EMAIL')
# Error timestamp notification window
ERROR_EMAIL_NOTIFICATION_WINDOW=300

# local data directory
DATA_DIR = "/usr/data/"
GENE_NAMES_FILE = DATA_DIR+'gene_names.txt'
GENOME_BUILDS_FILE = DATA_DIR+'genome_builds.txt'

# local executables
BLASTP_EXECUTABLE = "/usr/externals/blast/bin/blastp"
CLUSTALW_EXECUTABLE = "/usr/externals/clustalw/clustalw2"
HMMFETCH_EXECUTABLE = "/usr/externals/hmmer/binaries/hmmfetch"
HMMLOGO_EXECUTABLE = "/usr/externals/hmmer/binaries/hmmlogo"
HMMALIGN_EXECUTABLE = "/usr/externals/hmmer/binaries/hmmalign"
HMMEMIT_EXECUTABLE = "/usr/externals/hmmer/binaries/hmmemit"
HMMSTAT_EXECUTABLE = "/usr/externals/hmmer/binaries/hmmstat"

# Genome specific files
GENCODE_HG_ANNOTATION_FILE_GTF = DATA_DIR+"Gencode/gencode.v19.annotation.gtf"
GENCODE_HG_ANNOTATION_FILE_GFF3 = DATA_DIR+"Gencode/gencode.v19.annotation.gff3"
GENCODE_HG_TRANSCRIPTION_FILE = DATA_DIR+"Gencode/gencode.v19.pc_transcripts.fa"
GENCODE_HG_TRANSLATION_FILE = DATA_DIR+"Gencode/gencode.v19.pc_translations.fa"
GENCODE_REFSEQ_FILE = DATA_DIR+"Gencode/gencode.v19.metadata.RefSeq"
GENCODE_SWISSPROT_FILE = DATA_DIR+"Gencode/gencode.v19.metadata.SwissProt"
GENCODE_BASIC_FILE = DATA_DIR+"Gencode/ucsc.gencode.v19.wgEncodeGencodeBasic.txt"

# InterPro Files
INTERPROSCAN_DOCKER_IMAGE = "blaxterlab/interproscan:5.22-61.0"
INTERPROSCAN_DOCKER_VOLUME = 'metadom_interpro_temp'
INTERPROSCAN_EXECUTABLE = "interproscan.sh"
INTERPROSCAN_TEMP_DIR = '/usr/interpro_temp'
INTERPROSCAN_DOMAIN_DATABASES = 'Pfam'

# UNIPROT
UNIPROT_MAX_BLAST_RESULTS = 10
UNIPROT_DIR = DATA_DIR+"UniProt/"
UNIPROT_SPROT_CANONICAL_AND_ISOFORM = UNIPROT_DIR+"uniprot_sprot_canonical_and_varsplic.fasta"
UNIPROT_SPROT_ISOFORM = UNIPROT_DIR+"uniprot_sprot_varsplic.fasta"
UNIPROT_SPROT_CANONICAL = UNIPROT_DIR+"uniprot_sprot.fasta"
UNIPROT_SPROT_SPECIES_FILTER = "HUMAN"

# Meta-domain files
METADOMAIN_DIR = DATA_DIR+"metadomains/"
RECONSTRUCT_METADOMAINS = False
METADOMAIN_ALIGNMENT_FILE_NAME = 'metadomain_alignments' # Alignments are saved as: METADOMAIN_DIR+<Pfam_id>+'/'+METADOMAIN_ALIGNMENT_FILE_NAME
METADOMAIN_MAPPING_FILE_NAME = 'metadomain_mappings' # Mappings are saved as: METADOMAIN_DIR+<Pfam_id>+'/'+METADOMAIN_MAPPING_FILE_NAME
METADOMAIN_DETAILS_FILE_NAME = 'metadomain_details.json' # Details are saved as: METADOMAIN_DIR+<Pfam_id>+'/'+METADOMAIN_DETAILS_FILE_NAME
METADOMAIN_SNV_ANNOTATION_FILE_NAME = 'metadomain_snv_annotation' # Annotations are saved as: METADOMAIN_DIR+<Pfam_id>+'/'+METADOMAIN_SNV_ANNOTATION_FILE_NAME

# Pre-building visualization settings
PRE_BUILD_VISUALIZATION_DIR = DATA_DIR+"metadome_visualization/"
PRE_BUILD_VISUALIZATION_FILE_NAME = 'metadome_visualization.json' # Visualizations are saved as: PRE_BUILD_VISUALIZATION_DIR+<Transcript_id>+'/'+PRE_BUILD_VISUALIZATION_FILE_NAME
PRE_BUILD_VISUALIZATION_TASK_FILE_NAME = 'visualization_task'
PRE_BUILD_VISUALIZATION_ERROR_FILE_NAME = 'visualization_error'
METADOMAIN_CACHE_MAXSIZE = 128 # set at 128 as in 2027 build the meta-domains with more than 100 occurrences is 84+
PREBUILD_PRINT_LIST_EXAMPLES_CAP = 25  # how many example ids to print per category

# PFAM specific files
PFAM_DIR = DATA_DIR+"PFAM/Pfam30.0"
PFAM_ALIGNMENT_DIR = PFAM_DIR+"/alignment/"
PFAM_HMM_DAT = PFAM_DIR+"/Pfam-A.hmm.dat.gz"
PFAM_HMM = PFAM_DIR+"/Pfam-A.hmm"

# gnomAD specific files
GNOMAD_DIR = DATA_DIR + "gnomAD/"
GNOMAD_GRCH37_VCF_FILE = GNOMAD_DIR + "GRCh37" + "/" + "pass_gnomad.exomes.r2.0.2.sites.vcf.gz"
GNOMAD_GRCH38_VCF_FILE = GNOMAD_DIR + "GRCh38" + "/" + "gnomad.joint.v4.1.sites.exomes.vcf.gz"
GNOMAD_ACCEPTED_FILTERS = ['PASS']

# ClinVar specific files
CLINVAR_DIR = DATA_DIR + 'ClinVar/'
CLINVAR_GRCH37_VCF_FILE = CLINVAR_DIR + "GRCh37" + "/" + 'clinvar_20251006.vcf.gz'
CLINVAR_GRCH38_VCF_FILE = CLINVAR_DIR + "GRCh38" + "/" + 'clinvar_20251006.vcf.gz'
CLINVAR_CONSIDERED_CLINSIG = ['Pathogenic']