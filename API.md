# G.O.A.T API Documentation

Complete API reference for all HTTP endpoints in the G.O.A.T web application.

---

## Table of Contents

1. [Main Workflow Routes](#main-workflow-routes)
2. [Data Retrieval Routes](#data-retrieval-routes)
3. [Utility Routes](#utility-routes)
4. [Response Formats](#response-formats)
5. [Error Handling](#error-handling)

---

## Main Workflow Routes

### Home Page

```
GET /
POST /
```

**Description**: Main dashboard showing all datasets and workflow entry point

**POST Parameters**:
- `type` (string): Workflow type ('filter', 'bin', 'compare', 'setops', 'preprocess')
- Additional parameters depend on workflow type (see specific endpoints)

**Returns**: HTML page with flash messages

**Example**:
```bash
curl -X POST http://localhost:8000/ \
  -F "type=filter" \
  -F "dataset=original/sample1" \
  -F "save=filtered_output"
```

---

### Filter - Small RNA Species Selection

```
GET /filter
```

**Description**: Filter small RNA reads by length, nucleotide type, and custom rules

**Form Parameters** (POST to /):
- `type`: "filter"
- `dataset` (string): Path to dataset (e.g., "original/experiment1/sample1")
- `length_start` (int): Minimum sequence length (10-35)
- `length_end` (int): Maximum sequence length (10-35)
- `nucleotides` (list): Nucleotide types to include (['A'], ['T'], ['G'], ['C'], or combinations)
- `strand` (string): Strand selection ('s' for sense, 'a' for antisense)
- `hidden_rule_<group>_<id>` (string): Custom filtering rules (format: "operator;value")
- `save` (string): Output filename

**Rule Format**:
- Sequence-level: `>`, `<`, `=`, `≥`, `≤` followed by PPM threshold
- Gene-level: `Top` or `Bottom` followed by percentage

**Example**:
```python
{
    'type': 'filter',
    'dataset': 'original/exp1/sample',
    'length_start': '21',
    'length_end': '23',
    'nucleotides': ['G'],
    'strand': 's',
    'hidden_rule_0_1': '>;100',  # Sequences with PPM > 100
    'save': '21-23G_filtered'
}
```

**Returns**:
- Success: `"Saved as <filename>"`
- Error: `"Nothing found, file is not created!"`

---

### Bin - Gene-Level Aggregation

```
GET /bin
```

**Description**: Aggregate read counts to gene level with filtering

**Form Parameters** (POST to /):
- `type`: "bin"
- `dataset` (string): Path to filtered dataset
- `unique_mappers` (bool): Include unique mappers (checkbox)
- `multi_mappers` (bool): Include multi-mappers (checkbox)
- `gene_types` (list): Gene types to include (['protein_coding', 'miRNA', etc.])
- `gene_sets` (list): Predefined gene sets to filter by
- `hasSequence` (bool): Data contains sequence information
- `isTranscript` (bool): Data uses transcript IDs
- `save` (string): Output filename

**Example**:
```python
{
    'type': 'bin',
    'dataset': 'filtered/21-23G_filtered',
    'unique_mappers': 'on',
    'multi_mappers': 'on',
    'gene_types': ['protein_coding', 'miRNA'],
    'gene_sets': ['wago', 'csr'],
    'save': 'binned_output'
}
```

**Returns**: Statistics and filename

---

### Compare - Multi-Dataset Comparison

```
GET /compare
```

**Description**: Compare up to 4 datasets with complex rules

**Form Parameters** (POST to /):
- `type`: "compare"
- `groupA`, `groupB`, `groupC`, `groupD` (string): Paths to binned datasets
- `hidden_rule_<group>_<id>` (string): Comparison rules (format: "dataset;operator;value")
- `save` (string): Output gene list name

**Rule Format**: `<dataset_or_ratio>;<operator>;<threshold>`
- Datasets: A, B, C, D
- Ratios: A/B, A/C, B/C, etc.
- Operators: >, <, =, ≥, ≤
- Example: `A/B;>;2` (A/B ratio > 2)

**Example**:
```python
{
    'type': 'compare',
    'groupA': 'binned/mutant',
    'groupB': 'binned/wildtype',
    'hidden_rule_0_1': 'A/B;>;2',  # Mutant > 2x wildtype
    'hidden_rule_0_2': 'A;>;100',  # Mutant PPM > 100
    'save': 'enriched_in_mutant'
}
```

**Returns**: Number of genes and filename

---

### Set Operations - Venn Diagrams

```
GET /setops
```

**Description**: Perform set operations on 2-4 gene lists

**Form Parameters** (POST to /):
- `type`: "setops"
- `groupA`, `groupB`, `groupC`, `groupD` (string): Gene list names
- `hidden_operation` (string): Operation type
- `save` (string): Output gene list name

**Supported Operations**:
- **2-way**: `and` (A∩B), `or` (A∪B), `a-b` (A-B), `b-a` (B-A)
- **3-way**: `abc-and` (A∩B∩C), `abc-or` (A∪B∪C), `ab-c`, `ac-b`, `bc-a`, `a-bc`, `b-ac`, `c-ab`
- **4-way**: `abcd-and`, `abcd-or`, `abc-d`, `abd-c`, `acd-b`, `bcd-a`, `a-bcd`, `b-acd`, `c-abd`, `d-abc`

**Example**:
```python
{
    'type': 'setops',
    'groupA': 'wago_targets',
    'groupB': 'csr_targets',
    'groupC': 'oogenic',
    'hidden_operation': 'abc-and',  # Intersection of all 3
    'save': 'common_targets'
}
```

---

### Plot - Data Visualization

```
GET /plot
GET /plot/<plot_name>
```

**Description**: Generate interactive plots

**Available Plot Types** (from `templates/plot_templates/`):
- `basic_scatter` - Simple X vs Y scatter plot
- `fast_scatter_plot` - Optimized scatter plot
- Custom plots as added

**Returns**: HTML page with plot interface

---

### Preprocess - Raw Data Pipeline

```
GET /preprocess
POST /preprocess
```

**Description**: Process raw FASTQ files through alignment pipeline

**Form Parameters**:
- `path` (string): Directory containing FASTQ files
- `username` (string): User identifier
- `experiment_name` (string): Experiment name
- `filename_<N>` (string): FASTQ filename
- `alias_<N>` (string): Sample alias
- `use_counts_fa` (bool): Use counts FASTA workflow
- `adapter_seq` (string): Adapter sequence for trimming
- `min_seq_len` (int): Minimum sequence length after trimming
- `first_n` (int): Number of reads to process (0 for all)

**Example**:
```python
{
    'path': '/data/raw_fastq',
    'username': 'researcher1',
    'experiment_name': 'exp1',
    'filename_0': 'sample1.fastq.gz',
    'alias_0': 'wildtype',
    'adapter_seq': 'TGGAATTCTCGGGTGCCAAGG',
    'min_seq_len': '18',
    'first_n': '0'
}
```

**Process**:
1. Cutadapt: Adapter trimming
2. Bowtie2: Alignment to genome
3. SAM parsing: Quantification
4. Data splitting: By length/nucleotide

**Returns**: Empty string (processes in background)

---

## Data Retrieval Routes

### Preview Set Operation

```
GET /preview/<groupA>/<groupB>/<groupC>/<groupD>/<operation>
```

**Description**: Preview results of set operation without saving

**Parameters**:
- `groupA`, `groupB`, `groupC`, `groupD`: Gene list names (use `|` for empty)
- `operation`: Set operation type

**Returns**: JSON
```json
{
  "output": ["gene1", "gene2", ...],
  "list_a": ["gene1", ...],
  "list_b": ["gene2", ...],
  "list_c": [...],
  "list_d": [...]
}
```

**Example**:
```bash
curl http://localhost:8000/preview/wago/csr/|/|/and
```

---

### Get Data Pair for Plotting

```
GET /getDataPair/<fileA>/<fileB>/<folder>
```

**Description**: Retrieve gene expression data for two datasets

**Parameters**:
- `fileA`: Filename in folder
- `fileB`: Filename in folder
- `folder`: Folder path with underscores (e.g., `data_binned`)

**Returns**: JSON array
```json
[
  {
    "gene": "WBGene00001234",
    "type": "protein_coding",
    "x": 123.5,
    "y": 98.2
  },
  ...
]
```

**Example**:
```bash
curl http://localhost:8000/getDataPair/mutant.json/wildtype.json/data_binned
```

---

### Get Gene List

```
GET /getgenelist/<genelist>
```

**Description**: Retrieve genes in a gene list (converted to locus names)

**Returns**: JSON array of gene names
```json
["Y74C9A.6", "F56D12.5", ...]
```

---

### Get File List

```
GET /getlist?path=<path>&file_type=<type>
```

**Description**: List files in directory with specific extension

**Parameters**:
- `path`: Directory path
- `file_type`: File extension (`.fastq.gz`, `.fa`)

**Returns**: JSON array of filenames

---

## Utility Routes

### Discover Gene List Origins

```
GET /discover
```

**Description**: Explore which gene lists contain specific genes

**Returns**: HTML page with gene list analysis

---

### Venn Diagram Development

```
GET /venn
```

**Description**: Interactive Venn diagram development/testing page

---

## Response Formats

### Success Response (POST workflows)

**Flash Message**:
```
"Saved as <filename>, having <N> genes"
"<N> genes passing the filter. Saved as <filename>"
```

### Error Response

**Flash Message**:
```
"!No genes found after set operation, file is not saved."
"Nothing found, file is not created!"
"No genes passed the filter, data not saved."
```

### JSON Response (Data endpoints)

**Success**:
```json
{
  "data": [...],
  "metadata": {...}
}
```

**Error**:
```json
{
  "error": "Error message here"
}
```

---

## Error Handling

### HTTP Status Codes

- `200 OK`: Successful request
- `400 Bad Request`: Invalid parameters
- `404 Not Found`: Resource not found
- `500 Internal Server Error`: Server error

### Common Errors

1. **File Not Found**
   - Cause: Dataset path doesn't exist
   - Solution: Check file paths, use `/getlist` to verify

2. **Invalid Parameters**
   - Cause: Missing required fields or invalid values
   - Solution: Validate input types and ranges

3. **Empty Results**
   - Cause: Filtering criteria too strict
   - Solution: Relax filtering rules

4. **Permission Denied**
   - Cause: Cannot write to output directory
   - Solution: Check directory permissions

---

## Rate Limits

**Production Deployment**:
- Preprocessing: 5 requests per hour per IP
- General API: 100 requests per minute per IP

**Development**:
- No rate limits

---

## Authentication

**Current**: No authentication required

**Future**: May add API keys for production deployments

---

## Versioning

**Current Version**: v2.0

**URL Format**: No version in URL (will add `/api/v2/` prefix in future)

---

## Examples

### Complete Workflow Example

```python
import requests

BASE_URL = "http://localhost:8000"

# 1. Filter for 22G RNAs
filter_response = requests.post(f"{BASE_URL}/", data={
    'type': 'filter',
    'dataset': 'original/exp1/sample1',
    'length_start': '22',
    'length_end': '22',
    'nucleotides': ['G'],
    'strand': 's',
    'save': '22G_sense'
})

# 2. Bin to genes
bin_response = requests.post(f"{BASE_URL}/", data={
    'type': 'bin',
    'dataset': 'filtered/22G_sense',
    'unique_mappers': 'on',
    'gene_types': ['protein_coding'],
    'save': '22G_genes'
})

# 3. Compare mutant vs wildtype
compare_response = requests.post(f"{BASE_URL}/", data={
    'type': 'compare',
    'groupA': 'binned/22G_genes_mutant',
    'groupB': 'binned/22G_genes_wildtype',
    'hidden_rule_0_1': 'A/B;>;2',
    'save': 'enriched_in_mutant'
})

# 4. Get gene list
genes = requests.get(f"{BASE_URL}/getgenelist/enriched_in_mutant").json()
print(f"Found {len(genes)} enriched genes")
```

---

## Best Practices

1. **Always validate input parameters** before submitting
2. **Use descriptive filenames** for traceability
3. **Check file existence** with `/getlist` before processing
4. **Preview set operations** before saving large gene lists
5. **Monitor flash messages** for errors and warnings

---

## Support

For API issues or questions:
- GitHub Issues: https://github.com/ahmetrasit/goat/issues
- Documentation: See README.md

---

**Last Updated**: 2025-11-18
**Maintainer**: G.O.A.T Development Team
