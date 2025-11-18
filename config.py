"""
Configuration management for G.O.A.T (Gene-Oriented Analysis Tool)
"""
import os
from pathlib import Path

# Base directory
BASE_DIR = Path(__file__).resolve().parent

# Flask configuration
class Config:
    """Base configuration"""
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'dev-secret-key-change-in-production'
    SEND_FILE_MAX_AGE_DEFAULT = 0
    DEBUG = os.environ.get('FLASK_DEBUG', 'False').lower() == 'true'

# Data directories
DATA_DIR = BASE_DIR / 'data'
DATA_ORIGINAL = DATA_DIR / 'original'
DATA_FILTERED = DATA_DIR / 'filtered'
DATA_BINNED = DATA_DIR / 'binned'
DATA_GENELIST = DATA_DIR / 'genelist'

# Mapper files directory
MAPPERS_DIR = BASE_DIR / 'mappers'
GENOME_INDEX = MAPPERS_DIR / 'merged.ws274'
GENOME_FASTA = MAPPERS_DIR / 'genome'

# External tools paths (override via environment variables)
BOWTIE2_PATH = os.environ.get('BOWTIE2_PATH', '/usr/bin/bowtie2')
CUTADAPT_PATH = os.environ.get('CUTADAPT_PATH', 'cutadapt')
SAMTOOLS_PATH = os.environ.get('SAMTOOLS_PATH', 'samtools')
HISAT2_PATH = os.environ.get('HISAT2_PATH', 'hisat2')
UCSC_TOOLS_DIR = os.environ.get('UCSC_TOOLS_DIR', '/usr/local/bin')

# Processing parameters
PROCESSING_THREADS = int(os.environ.get('PROCESSING_THREADS', os.cpu_count() // 2))
WORKER_POOL_SIZE = int(os.environ.get('WORKER_POOL_SIZE', os.cpu_count() // 2))

# Genome information
GENOME_VERSION = 'ws274'
SPECIES = 'Caenorhabditis elegans'

# Small RNA parameters
DEFAULT_MIN_LENGTH = 10
DEFAULT_MAX_LENGTH = 35
DEFAULT_ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'  # Illumina small RNA adapter
DEFAULT_MIN_ADAPTER_LENGTH = 18

# ID conversion mappers
ID_MAPPERS = {
    'alias2name': 'alias2name.json',
    'name2gene': 'name2gene.json',
    'gene2name': 'gene2name.json',
    'alias2type': 'alias2type.json',
    'name2type': 'name2type.json',
    'gene2type': 'gene2type.json',
}

# Gene types (WormBase biotypes)
GENE_TYPES = [
    'protein_coding',
    'miRNA',
    'ncRNA',
    'pseudogene',
    'transposon',
    'rRNA',
    'lincRNA',
    'piRNA',
    'tRNA',
    'snoRNA',
    'asRNA',
    'snRNA',
    'scRNA',
]

# Logging configuration
LOG_LEVEL = os.environ.get('LOG_LEVEL', 'INFO')
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

def ensure_directories():
    """Create necessary directories if they don't exist"""
    directories = [
        DATA_DIR,
        DATA_ORIGINAL,
        DATA_FILTERED,
        DATA_BINNED,
        DATA_GENELIST,
        MAPPERS_DIR,
    ]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)

if __name__ == '__main__':
    # Validate configuration
    print("G.O.A.T Configuration")
    print("=" * 50)
    print(f"Base Directory: {BASE_DIR}")
    print(f"Data Directory: {DATA_DIR}")
    print(f"Mappers Directory: {MAPPERS_DIR}")
    print(f"Bowtie2: {BOWTIE2_PATH}")
    print(f"Processing Threads: {PROCESSING_THREADS}")
    print(f"Genome Version: {GENOME_VERSION}")
    print("=" * 50)
    ensure_directories()
    print("✓ All directories verified/created")
