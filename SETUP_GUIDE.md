# G.O.A.T Installation and Setup Guide

Complete step-by-step guide for installing and running G.O.A.T (Gene-Oriented Analysis Tool) for C. elegans small RNA-seq analysis.

---

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Quick Start](#quick-start)
3. [Detailed Installation](#detailed-installation)
4. [Configuration](#configuration)
5. [Running the Application](#running-the-application)
6. [Testing](#testing)
7. [Troubleshooting](#troubleshooting)
8. [Production Deployment](#production-deployment)

---

## System Requirements

### Hardware

- **CPU**: Multi-core processor (4+ cores recommended)
- **RAM**: 8 GB minimum, 16+ GB recommended
- **Storage**: 100+ GB for genome indices and data
- **OS**: Linux/Unix (Ubuntu 20.04+ or similar)

### Software Prerequisites

- **Python**: 3.8 or higher
- **Git**: For version control
- **Build tools**: GCC, make
- **Internet connection**: For downloading dependencies and genome data

---

## Quick Start

For experienced users who want to get started immediately:

```bash
# Clone and setup
git clone https://github.com/ahmetrasit/goat.git
cd goat
git checkout claude/code-review-implementation-01HPTQ4ept53RKtzNMty9mdn

# Install system dependencies
sudo apt-get update
sudo apt-get install -y bowtie2 cutadapt samtools hisat2 bedtools

# Setup Python environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Configure
cp .env.example .env
# Edit .env with your paths

# Run
python app.py
# Access at http://localhost:8000
```

---

## Detailed Installation

### Step 1: Install System Dependencies

#### Ubuntu/Debian

```bash
# Update package list
sudo apt-get update

# Install bioinformatics tools
sudo apt-get install -y \
    bowtie2 \
    cutadapt \
    samtools \
    hisat2 \
    bedtools \
    build-essential \
    python3-dev \
    python3-pip \
    python3-venv

# Verify installations
bowtie2 --version
cutadapt --version
samtools --version
```

#### UCSC Tools

```bash
# Create directory for UCSC tools
sudo mkdir -p /usr/local/bin/ucsc-tools
cd /usr/local/bin/ucsc-tools

# Download tools (Linux x86_64)
sudo wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
sudo wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig

# Make executable
sudo chmod +x *

# Add to PATH (add to ~/.bashrc for permanent)
export PATH=$PATH:/usr/local/bin/ucsc-tools
```

---

### Step 2: Clone Repository

```bash
# Clone the repository
git clone https://github.com/ahmetrasit/goat.git
cd goat

# Checkout the implementation branch
git checkout claude/code-review-implementation-01HPTQ4ept53RKtzNMty9mdn

# Verify you're on the correct branch
git branch
```

---

### Step 3: Setup Python Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install Python dependencies
pip install -r requirements.txt

# Verify installation
python -c "import flask; print(f'Flask {flask.__version__} installed')"
```

---

### Step 4: Download Genome Data

#### Option A: C. elegans WS274 (Recommended)

```bash
# Create mappers directory
mkdir -p mappers
cd mappers

# Download genome from WormBase
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS274.genomic.fa.gz

# Extract
gunzip c_elegans.PRJNA13758.WS274.genomic.fa.gz

# Build Bowtie2 index (takes 10-20 minutes)
bowtie2-build c_elegans.PRJNA13758.WS274.genomic.fa merged.ws274

# Build HISAT2 index for JBrowse (optional)
hisat2-build c_elegans.PRJNA13758.WS274.genomic.fa genome

# Verify index files exist
ls -lh merged.ws274.*.bt2
```

#### Gene Annotation Files

The repository includes pre-built annotation files in `mappers/`:
- `gene2name.json`
- `name2type.json`
- `alias2gene.json`
- etc.

If these are missing, you can regenerate them from WormBase GFF files (see documentation).

---

### Step 5: Setup Data Directories

```bash
# Return to project root
cd ..

# Create data directories
mkdir -p data/original
mkdir -p data/filtered
mkdir -p data/binned
mkdir -p data/genelist

# Set permissions
chmod 755 data
chmod -R 755 data/*

# Verify structure
tree -L 2 data/
```

---

## Configuration

### Environment Variables

```bash
# Copy template
cp .env.example .env

# Edit configuration
nano .env
```

**Required Settings**:

```bash
# Security (IMPORTANT: Generate a strong random key for production)
SECRET_KEY=change-this-to-a-random-secret-key

# External Tools Paths
BOWTIE2_PATH=/usr/bin/bowtie2
CUTADAPT_PATH=/usr/bin/cutadapt
SAMTOOLS_PATH=/usr/bin/samtools
HISAT2_PATH=/usr/bin/hisat2
UCSC_TOOLS_DIR=/usr/local/bin/ucsc-tools

# Processing Configuration
PROCESSING_THREADS=8  # Adjust based on your CPU cores
WORKER_POOL_SIZE=4

# Logging
LOG_LEVEL=INFO
```

**Generate Secret Key**:

```python
# In Python shell
import secrets
print(secrets.token_hex(32))
# Copy output to SECRET_KEY in .env
```

---

## Running the Application

### Development Mode

```bash
# Activate virtual environment (if not already)
source venv/bin/activate

# Run development server
python app.py

# Application will be available at:
# http://localhost:8000
```

**Development Features**:
- Auto-reload on code changes
- Debug mode enabled
- Detailed error pages

---

### Production Mode

```bash
# Activate virtual environment
source venv/bin/activate

# Run with Gunicorn
gunicorn -c gunicorn_config.py app:app

# Or with custom workers
gunicorn -w 4 -b 0.0.0.0:8000 app:app
```

**Production Features**:
- 20 worker processes (configurable)
- Auto-restart on failures
- Production-grade logging
- Security headers enabled

---

### Running as System Service

Create systemd service file:

```bash
sudo nano /etc/systemd/system/goat.service
```

```ini
[Unit]
Description=G.O.A.T Small RNA Analysis Tool
After=network.target

[Service]
Type=notify
User=www-data
Group=www-data
WorkingDirectory=/path/to/goat
Environment="PATH=/path/to/goat/venv/bin"
ExecStart=/path/to/goat/venv/bin/gunicorn -c gunicorn_config.py app:app
Restart=always

[Install]
WantedBy=multi-user.target
```

```bash
# Reload systemd
sudo systemctl daemon-reload

# Start service
sudo systemctl start goat

# Enable on boot
sudo systemctl enable goat

# Check status
sudo systemctl status goat
```

---

## Testing

### Run Test Suite

```bash
# Activate virtual environment
source venv/bin/activate

# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov --cov-report=html

# View coverage report
open htmlcov/index.html  # Or xdg-open on Linux
```

### Test Specific Modules

```bash
# Test set operations
pytest tests/test_process.py -v

# Test security functions
pytest tests/test_utils.py::TestPathValidation -v

# Test SAM parsing
pytest tests/test_transcript_analysis.py -v
```

### Expected Output

```
============================= test session starts ==============================
collected 25 items

tests/test_process.py::TestSetOperations::test_applyOperation_2way_and PASSED
tests/test_process.py::TestSetOperations::test_applyOperation_3way_and PASSED
tests/test_process.py::TestSetOperations::test_applyOperation_4way_and PASSED
...

---------- coverage: platform linux, python 3.8.10 -----------
Name                          Stmts   Miss  Cover
-------------------------------------------------
process.py                      250     25    90%
utils.py                        180     15    92%
TranscriptAnalysis.py           300     60    80%
-------------------------------------------------
TOTAL                           730    100    86%
```

---

## Troubleshooting

### Common Issues

#### 1. Port Already in Use

**Error**: `Address already in use`

**Solution**:
```bash
# Find process using port 8000
sudo lsof -i :8000

# Kill the process
sudo kill -9 <PID>

# Or use a different port
python app.py --port 8080
```

#### 2. Permission Denied

**Error**: `Permission denied: 'data/binned/output.json'`

**Solution**:
```bash
# Fix data directory permissions
chmod -R 755 data/
chown -R $USER:$USER data/
```

#### 3. Bowtie2 Index Not Found

**Error**: `FileNotFoundError: mappers/merged.ws274.1.bt2`

**Solution**:
```bash
# Rebuild index
cd mappers
bowtie2-build c_elegans.PRJNA13758.WS274.genomic.fa merged.ws274
cd ..
```

#### 4. Import Errors

**Error**: `ModuleNotFoundError: No module named 'flask'`

**Solution**:
```bash
# Ensure virtual environment is activated
source venv/bin/activate

# Reinstall dependencies
pip install -r requirements.txt
```

#### 5. UCSC Tools Not Found

**Error**: `FileNotFoundError: wigToBigWig`

**Solution**:
```bash
# Set path in .env
UCSC_TOOLS_DIR=/path/to/ucsc-tools

# Or add to system PATH
export PATH=$PATH:/usr/local/bin/ucsc-tools
```

---

## Production Deployment

### Security Checklist

- [ ] Generate strong `SECRET_KEY`
- [ ] Set `FLASK_DEBUG=False` in production
- [ ] Configure HTTPS (use nginx/Apache as reverse proxy)
- [ ] Set up firewall rules
- [ ] Configure rate limiting with Redis
- [ ] Set proper file permissions (644 for files, 755 for dirs)
- [ ] Disable directory listing
- [ ] Set up log rotation
- [ ] Configure backup strategy

### Nginx Reverse Proxy

```nginx
server {
    listen 80;
    server_name your-domain.com;

    # Redirect to HTTPS
    return 301 https://$server_name$request_uri;
}

server {
    listen 443 ssl http2;
    server_name your-domain.com;

    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location /static {
        alias /path/to/goat/static;
        expires 30d;
    }
}
```

### Redis for Rate Limiting (Optional)

```bash
# Install Redis
sudo apt-get install redis-server

# Start Redis
sudo systemctl start redis
sudo systemctl enable redis

# Update .env
RATE_LIMIT_STORAGE_URI=redis://localhost:6379

# Uncomment in requirements.txt
# redis>=4.6.0
```

---

## Verification

### Health Checks

```bash
# Check application is running
curl http://localhost:8000/

# Check API endpoints
curl http://localhost:8000/filter
curl http://localhost:8000/bin

# Check file listing works
curl "http://localhost:8000/getlist?path=data/genelist&file_type=.json"
```

### Log Files

```bash
# Check application logs
tail -f logs/goat.log

# Check gunicorn logs
tail -f logs/gunicorn_error.log
```

---

## Next Steps

1. **Upload Sample Data**: Place FASTQ files in `data/original/`
2. **Process First Sample**: Use web interface to preprocess → filter → bin → compare
3. **Explore API**: See `API.md` for programmatic access
4. **Read Documentation**: Check `README.md` for detailed usage
5. **Run Tests**: Ensure everything works with `pytest tests/`

---

## Support

### Documentation
- **README.md**: User guide and features
- **API.md**: Complete API reference
- **code-review-recommendations.md**: Implementation details

### Getting Help
- **Issues**: https://github.com/ahmetrasit/goat/issues
- **Discussions**: GitHub Discussions
- **Email**: [Contact info]

---

## Quick Reference

### Start Application
```bash
source venv/bin/activate
python app.py  # Development
gunicorn -c gunicorn_config.py app:app  # Production
```

### Run Tests
```bash
pytest tests/ -v --cov
```

### Update Dependencies
```bash
pip install --upgrade -r requirements.txt
```

### Check Logs
```bash
tail -f logs/goat.log
```

### Backup Data
```bash
tar -czf goat_backup_$(date +%Y%m%d).tar.gz data/
```

---

**Version**: 2.0
**Last Updated**: 2025-11-18
**Maintained by**: G.O.A.T Development Team

🎉 **You're all set! Start analyzing small RNAs!** 🎉
