# Flask settings
from metadome.flask_app_credentials import SECRET_KEY_CRED
DEBUG = False
SECRET_KEY = SECRET_KEY_CRED

# FLask-SQLAchemy settings
# from metadome.postgres_credentials import POSTGRES_USER, POSTGRES_PASSWORD
import os
SQLALCHEMY_RECORD_QUERIES = DEBUG # should be false when not debug
SQLALCHEMY_TRACK_MODIFICATIONS = False
SQLALCHEMY_DATABASE_URI = "postgresql://"+os.environ["POSTGRES_USER"]+":"+os.environ["POSTGRES_PASSWORD"]+"@"+os.environ["POSTGRES_HOST"]+"/"+os.environ["POSTGRES_DB"]
SQLALCHEMY_ECHO = True
SQLALCHEMY_POOL_TIMEOUT = 10

# Flask-Celery settings
CELERY_BROKER_URL='amqp://guest@metadome-rabbitmq-1'
CELERY_RESULT_BACKEND='redis://metadome-redis-1/0'
CELERY_TRACK_STARTED = True
CELERY_TASK_SERIALIZER = 'pickle'
CELERY_RESULT_SERIALIZER='pickle'
CELERY_ACCEPT_CONTENT = ['pickle']

# Visualiation specific settings
ALLELE_FREQUENCY_CUTOFF = 0.0
SLIDING_WINDOW_SIZE = 10

# Debug toolbar
DEBUG_TB_ENABLED = DEBUG

# E-mail
MAIL_SERVER = None # add your smtp server here if needed
DEFAULT_RECIPIENT = None # where should the support emails be sent to

# local data directory
DATA_DIR = "/usr/data/"
GENE_NAMES_FILE = DATA_DIR+'gene_names.txt'

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

# Pre-build visualization files
PRE_BUILD_VISUALIZATION_DIR = DATA_DIR+"metadome_visualization/"
PRE_BUILD_VISUALIZATION_FILE_NAME = 'metadome_visualization.json' # Visualizations are saved as: PRE_BUILD_VISUALIZATION_DIR+<Transcript_id>+'/'+PRE_BUILD_VISUALIZATION_FILE_NAME
PRE_BUILD_VISUALIZATION_TASK_FILE_NAME = 'visualization_task'
PRE_BUILD_VISUALIZATION_ERROR_FILE_NAME = 'visualization_error'

# PFAM specific files
PFAM_DIR = DATA_DIR+"PFAM/Pfam30.0"
PFAM_ALIGNMENT_DIR = PFAM_DIR+"/alignment/"
PFAM_HMM_DAT = PFAM_DIR+"/Pfam-A.hmm.dat.gz"
PFAM_HMM = PFAM_DIR+"/Pfam-A.hmm"

# gnomAD specific files
GNOMAD_DIR = DATA_DIR + "gnoMAD/"
GNOMAD_VCF_FILE = GNOMAD_DIR + "pass_gnomad.exomes.r2.0.2.sites.vcf.gz"
GNOMAD_ACCEPTED_FILTERS = ['PASS']

# ClinVar specific files
CLINVAR_DIR = DATA_DIR + 'ClinVar/GRCh37/'
CLINVAR_ORIGINAL_VCF_FILE = CLINVAR_DIR + 'clinvar.vcf.gz'
CLINVAR_VCF_FILE = CLINVAR_ORIGINAL_VCF_FILE
CLINVAR_CONSIDERED_CLINSIG = ['Pathogenic']
