"""
Constants and configuration values for G.O.A.T
Extracted from magic numbers throughout the codebase for better maintainability
"""

# Small RNA Sequence Parameters
MIN_SEQUENCE_LENGTH = 10  # Minimum valid small RNA length (nt)
MAX_SEQUENCE_LENGTH = 35  # Maximum valid small RNA length (nt)
MIN_VALID_SEQ_LENGTH_FOR_MAPPING = 10  # Sequences must be at least this long to be considered for mapping

# Small RNA Species Definitions
# Based on C. elegans small RNA biology
SMALL_RNA_SPECIES = {
    '21U': {'length': 21, 'nucleotide': 'U', 'description': 'piRNAs (21U RNAs)'},
    '22G': {'length': 22, 'nucleotide': 'G', 'description': 'Secondary siRNAs (WAGO pathway)'},
    '26G': {'length': 26, 'nucleotide': 'G', 'description': 'Primary siRNAs (CSR-1 pathway)'},
    'miRNA': {'length': 22, 'nucleotide': None, 'description': 'MicroRNAs'},
}

# Normalization
PPM_FACTOR = 1_000_000  # Parts per million normalization factor

# Processing Parameters
DEFAULT_CPU_FRACTION = 0.5  # Use half of available CPUs by default
MIN_WORKERS = 1
MAX_WORKERS = 32

# File Naming
MAX_FILENAME_RENAME_ATTEMPTS = 1000  # Prevent infinite loops in file naming
FILENAME_COUNTER_START = 1

# Gene Level Filtering
# Used in furtherFilter for top/bottom percentage selection
PERCENTAGE_SCALE = 100
MIN_GENE_RULE_PERCENTAGE = 0
MAX_GENE_RULE_PERCENTAGE = 100

# Multi-mapper Categories
MULTI_MAPPER_CATEGORIES = ['all', 'some', 'none']
UNIQUE_MAPPER_KEY = 'u'
MULTI_MAPPER_KEY = 'm'

# Strand Notation
STRAND_SENSE = 's'
STRAND_ANTISENSE = 'a'
STRAND_FULL_NAMES = {
    STRAND_SENSE: 'sense',
    STRAND_ANTISENSE: 'antisense',
    's': 'sense',
    'a': 'asense',  # Abbreviated form used in filenames
}

# SAM File Parsing
SAM_FLAG_FORWARD = 0
SAM_FLAG_REVERSE = 16
SAM_FLAG_MODULO = 256

# Gene Types (WormBase biotypes)
GENE_TYPES_DISPLAY = [
    'PC',           # Protein coding
    'MIRNA',        # MicroRNA
    'NCRNA',        # Non-coding RNA
    'PSEUODOGENE',  # Pseudogene (note: typo in original, keeping for consistency)
    'TRANSPOSON',   # Transposon
    'RRNA',         # Ribosomal RNA
    'LINCRNA',      # Long intergenic non-coding RNA
    'PIRNA',        # PIWI-interacting RNA
    'TRNA',         # Transfer RNA
    'SNORNA',       # Small nucleolar RNA
    'ASRNA',        # Antisense RNA
    'SNRNA',        # Small nuclear RNA
    'SCRNA',        # Small cytoplasmic RNA
]

# Structural RNA types to exclude in some analyses
STRUCTURAL_RNA_TYPES = {'ASRNA', 'LINCRNA', 'RRNA', 'SCRNA', 'SNORNA', 'SNRNA'}

# Data Folders
FOLDER_TYPES = ['original', 'filtered', 'binned', 'genelist']

# Precision for PPM rounding
PPM_DECIMAL_PLACES = 4

# Default adapter sequence for Illumina small RNA
DEFAULT_ILLUMINA_ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'

# Pseudocount for fold change calculations
DEFAULT_PSEUDOCOUNT = 1.0

# Comparison operators
COMPARISON_OPERATORS = ['<', '≤', '≥', '>', '=']
OPERATOR_SYMBOLS = {
    '<': lambda a, b: a < b,
    '≤': lambda a, b: a <= b,
    '≥': lambda a, b: a >= b,
    '>': lambda a, b: a > b,
    '=': lambda a, b: a == b,
}

# Set Operations
SET_OPERATIONS_2WAY = ['and', 'or', 'a-b', 'b-a']
SET_OPERATIONS_3WAY = ['abc-and', 'abc-or', 'ab-c', 'ac-b', 'bc-a', 'a-bc', 'b-ac', 'c-ab']
SET_OPERATIONS_4WAY = [
    'abcd-and', 'abcd-or',
    'abc-d', 'abd-c', 'acd-b', 'bcd-a',
    'a-bcd', 'b-acd', 'c-abd', 'd-abc'
]

# Validation
MAX_DATASET_NAME_LENGTH = 255
VALID_NUCLEOTIDES = {'A', 'T', 'G', 'C', 'U'}
VALID_FILE_EXTENSIONS = {'.json', '.fa', '.fasta', '.fastq', '.gz'}

# HTTP Status Codes
HTTP_OK = 200
HTTP_BAD_REQUEST = 400
HTTP_NOT_FOUND = 404
HTTP_INTERNAL_ERROR = 500

# Rate Limiting (requests per time period)
RATE_LIMIT_PREPROCESS = "5 per hour"  # Preprocessing is resource-intensive
RATE_LIMIT_GENERAL = "100 per minute"  # General API calls
RATE_LIMIT_STORAGE_URI = "redis://localhost:6379"  # For production

# Security
MIN_PASSWORD_LENGTH = 8  # If authentication is added
ALLOWED_UPLOAD_EXTENSIONS = {'.fastq.gz', '.fa', '.fasta'}
MAX_UPLOAD_SIZE_MB = 1000  # Maximum upload size in MB

# Logging
LOG_FORMAT_DETAILED = '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
LOG_FORMAT_SIMPLE = '%(asctime)s - %(levelname)s - %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'

# Performance
LARGE_FILE_THRESHOLD_MB = 100  # Files larger than this should use streaming
CHUNK_SIZE_BYTES = 8192  # For file streaming

# Genome Information
GENOME_ORGANISMS = {
    'ce': 'Caenorhabditis elegans',
    'cb': 'Caenorhabditis briggsae',
}
DEFAULT_ORGANISM = 'ce'
DEFAULT_GENOME_VERSION = 'WS274'

# ID Conversion Confidence
ID_CONVERSION_HIGH_CONFIDENCE = 90  # Percentage threshold for confident ID type detection
ID_CONVERSION_MEDIUM_CONFIDENCE = 50
ID_CONVERSION_LOW_CONFIDENCE = 25

# Plot Parameters
DEFAULT_PLOT_WIDTH = 800
DEFAULT_PLOT_HEIGHT = 600
SCATTER_POINT_SIZE = 3
HEATMAP_COLOR_SCHEME = 'RdYlBu'

# Error Messages
ERROR_FILE_NOT_FOUND = "File not found: {}"
ERROR_INVALID_JSON = "Invalid JSON format: {}"
ERROR_INVALID_PARAMETERS = "Invalid parameters: {}"
ERROR_NO_DATA = "No data found matching the criteria"
ERROR_PERMISSION_DENIED = "Permission denied: {}"

# Success Messages
SUCCESS_FILE_SAVED = "File saved successfully: {}"
SUCCESS_PROCESSING_COMPLETE = "Processing completed: {} genes"
SUCCESS_DATA_LOADED = "Data loaded: {} items"
