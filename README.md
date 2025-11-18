# G.O.A.T - Gene-Oriented Analysis Tool

A web-based bioinformatics pipeline for analyzing small RNA sequencing data in *Caenorhabditis elegans*.

## Features

- **Preprocessing Pipeline**: Adapter trimming, alignment, and quantification
- **Flexible Filtering**: Filter by length, nucleotide type, strand, and custom rules
- **Gene-level Binning**: Aggregate read counts to gene level
- **Multi-dataset Comparison**: Compare up to 4 datasets with complex rules
- **Set Operations**: Perform 2-4 way Venn diagrams on gene lists
- **Interactive Visualizations**: Scatter plots, gene type comparisons, heatmaps
- **JBrowse Integration**: Generate browser tracks for visualization

## Installation

### Prerequisites

**System Requirements:**
- Linux/Unix system
- Python 3.8+
- 8+ GB RAM recommended
- 100+ GB storage for genome indices and data

**External Tools:**
```bash
# Install bioinformatics tools (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y \
    bowtie2 \
    cutadapt \
    samtools \
    hisat2 \
    bedtools

# UCSC tools (download binaries)
# Visit: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
```

### Python Setup

```bash
# Clone repository
git clone <repository-url>
cd goat

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt

# Copy and configure environment
cp .env.example .env
# Edit .env with your paths
```

### Genome Data Setup

1. **Download C. elegans genome** (WS274):
```bash
mkdir -p mappers
cd mappers

# Download from WormBase
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS274.genomic.fa.gz

# Build Bowtie2 index
bowtie2-build c_elegans.PRJNA13758.WS274.genomic.fa.gz merged.ws274
```

2. **Prepare gene annotation mappers** (provided in `mappers/` directory):
- `gene2name.json` - Gene ID to gene name
- `name2type.json` - Gene name to biotype
- `alias2gene.json` - Alias to gene ID
- (See documentation for complete list)

## Usage

### Start the Server

```bash
# Development
python app.py

# Production (with gunicorn)
gunicorn -c gunicorn_config.py app:app
```

Access the web interface at `http://localhost:8000`

### Workflow

1. **Preprocess** - Upload FASTQ files, trim adapters, align to genome
2. **Filter** - Select small RNA species (21G, 22G, 26G, piRNAs, etc.)
3. **Bin** - Aggregate to gene level, filter by gene type/sets
4. **Compare** - Find differentially expressed genes
5. **Set Operations** - Combine gene lists with Venn diagrams
6. **Plot** - Visualize results

### Example Analysis

```python
# Example: Analyze 21G RNAs in sense orientation
# 1. Preprocess raw data
# 2. Filter:
#    - Length: 21-21 nt
#    - Nucleotide: G
#    - Strand: sense
# 3. Bin to genes (unique + multi-mappers)
# 4. Compare mutant vs wildtype (threshold: A/B > 2)
# 5. Intersect with WAGO targets (Set Operations)
```

## Configuration

### Environment Variables

See `.env.example` for all options:

```bash
SECRET_KEY=your-secret-key-here
BOWTIE2_PATH=/usr/bin/bowtie2
PROCESSING_THREADS=8
UCSC_TOOLS_DIR=/usr/local/bin
```

### File Structure

```
goat/
├── app.py                 # Flask application
├── process.py             # Data processing engine
├── preprocess.py          # Preprocessing pipeline
├── TranscriptAnalysis.py  # SAM file analysis
├── plot.py                # Plotting utilities
├── config.py              # Configuration
├── templates/             # HTML templates
├── static/                # CSS, JS, images
├── data/                  # Data storage (gitignored)
│   ├── original/          # Raw preprocessed data
│   ├── filtered/          # Filtered datasets
│   ├── binned/            # Gene-level counts
│   └── genelist/          # Gene lists
└── mappers/               # Gene annotations & indices
```

## Scientific Background

G.O.A.T is designed for analyzing small RNA pathways in *C. elegans*, particularly:

- **21U RNAs (piRNAs)**: Germline-specific, ~21nt, U-rich
- **22G RNAs**: Secondary siRNAs, WAGO pathway
- **26G RNAs**: Primary siRNAs, CSR-1 pathway
- **miRNAs**: MicroRNAs, ~22nt

### Multi-mapper Handling

G.O.A.T categorizes multi-mapping reads:
- **'all'**: Maps to genes that all have unique mappers (proportional allocation)
- **'some'**: Mixed situation (weighted allocation)
- **'none'**: Maps only to genes without unique mappers (equal distribution)

### Normalization

- **PPM (Parts Per Million)**: Reads normalized to total mapped reads
- Optional: miRNA/tRNA normalization for specific analyses

## Troubleshooting

**Issue**: `FileNotFoundError: mappers/merged.ws274.1.bt2`
**Solution**: Run `bowtie2-build` to create genome index

**Issue**: Import error for `config`
**Solution**: Ensure `config.py` is in the same directory as `app.py`

**Issue**: Permission denied for external tools
**Solution**: Verify tool paths in `.env` and check executable permissions

## Citation

If you use G.O.A.T in your research, please cite:

```
[Citation information to be added]
```

## License

[License information to be added]

## Contact

[Contact information to be added]

---

**Version**: 2.0 (November 2025)
**Genome**: *C. elegans* WS274
**Maintained by**: [Maintainer information]
