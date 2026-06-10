# 16S-rRNA-Pipeline: K-mer Based Taxonomic Profiling

16S-rRNA-Pipeline is a Python package for high-throughput 16S rRNA amplicon sequencing analysis. It uses a K-mer voting algorithm and SQLite-based reference database to generate fast and reproducible genus-level abundance profiles from raw FASTQ files.

## Installation

The package can be installed via GitHub using pip. 

``` bash
pip install git+<https://github.com/XiaoqingHan/16S-rRNA-Pipeline.git>
```

For development mode:

``` bash
pip install -e .
```

## Quick Start

After installation, all tools are available as command-line interfaces. Each tool supports `-h` for help.

### 1. Build reference database

Convert a reference FASTA into a searchable K-mer database.

``` bash
16s-build-db -i ref/HOMD.fasta -d ref/ref_database.db -k 25
```

### 2. Preprocessing & Quality Control

Trim primers, perform quality filtering, and merge paired-end reads.

Paired-end:

``` bash
16s-preprocess -i test_data/paired -o output/cleaned_data -m paired --trim_primer --p1 <FORWARD_PRIMER> --p2 <REVERSE_PRIMER>
```

Single-end:
``` bash
16s-preprocess -i test_data/single -o output/cleaned_data -m single --skip_trim
```

### 3. Chimera Removal

Remove potential chimeric sequences using reference-based filtering.

``` bash
16s-rm-chimera -i output/cleaned_data -o output/no_chimera --db ref/ref_database.db -k 25
```

### 4. Taxonomic Profiling

Generate genus-level abundance profiles using K-mer voting against the reference database.

``` bash
16s-profile -i output/no_chimera -o output/abundance_out --db ref/ref_database.db -k 25
```

### 5. Summary & Visualization

Merge individual sample profiles into a single abundance matrix and generate composition plots.

``` bash
16s-summarize -i output/abundance_out -o output/final_report --top_n 20
```

## Data & Test Datasets

### Reference database
HOMD (Human Oral Microbiome Database), version HOMD_16S_rRNA_RefSeq_V16.02_full.

### Test datasets
#### Single-end dataset
* Accession: PRJEB86033  
* Region: V3–V4  
* Note: Primer sequences were removed prior to analysis

#### Paired-end dataset
* Accession: PRJNA555320  
* Region: V4  
* Forward primer: GTGCCAGCMGCCGCGGTAA    
* Reverse primer: GGACTACHVGGGTWTCTAAT  

### Input format
- Single-end: *.fastq.gz
- Paired-end: *_1.fastq.gz and *_2.fastq.gz (or *_R1.fastq.gz and *_R2.fastq.gz)

Each sample should be in its own directory or follow a consistent naming convention.

### Outputs

- `all_samples_genus_matrix.csv`: Genus abundance matrix across all samples.
- `composition.png`: A high-resolution stacked bar plot showing the top $N$ genera.
- `*.log`: Detailed processing logs for reproducibility.

## Validation

The pipeline has been validated using oral microbiome 16S rRNA datasets.

Validation includes:
* Single-end and paired-end datasets
* Multiple sample subsets
* Execution in clean Conda environments
* End-to-end pipeline validation from raw FASTQ files to abundance matrix and visualization outputs
  
## Requirements

- **Python**: 3.8+
- **Core Dependencies**: Pandas, Numpy, Matplotlib, Seaborn, Scipy.
- **Hardware**: Multi-threading is supported; 16GB+ RAM is recommended for large K-mer databases.

## Tips

- Use consistent FASTQ naming to avoid silent file skipping
- Increase k-mer size improves specificity but increases runtime
- Always use the same reference database version for comparisons

## Contact

#### Author: Xiaoqing Han

#### Email: [xhan723\@hotmail.com](mailto:xhan723@hotmail.com)
