# 16S-rRNA-Pipeline: K-mer Based Taxonomic Profiling

16S-rRNA-Pipeline is a comprehensive Python package designed for high-throughput 16S rRNA amplicon sequencing analysis. By leveraging a custom K-mer voting algorithm and SQLite-backed reference databases, it provides a fast and reproducible workflow from raw FASTQ files to genus-level abundance matrices.

## Installation

This package can be installed directly from GitHub using pip. We recommend using a virtual environment.

``` bash
pip install git+<https://github.com/XiaoqingHan/16S-rRNA-Pipeline.git>
```

*Note: For developers who wish to modify the source code, use `pip install -e .` after cloning the repository locally.*

## Validation

The current version has been validated on oral microbiome 16S rRNA sequencing datasets used in this study.

Validation tests included:

* Paired-end sequencing data
* Single-end sequencing data
* Multiple randomly selected samples
* Installation and execution in a clean Conda environment using the GitHub installation method

The pipeline has been verified from installation through generation of abundance profiles and summary reports.

## Quick Start Guide

Once installed, the pipeline provides several command-line tools. You do not need to run Python scripts manually. Each tool supports the `-h` option to display detailed usage instructions.

### 1. Prepare Reference Database

Convert your reference FASTA (e.g., HOMD) into a searchable K-mer database.

``` bash
16s-build-db -i HOMD.fasta -d ref_database.db -k 25
```

### 2. Preprocessing & Quality Control

Trim primers, filter low-quality reads, and merge paired-ends.

Examples for paired-end and single-end data, with or without primer trimming.

``` bash
16s-preprocess -i ./raw_data -o ./cleaned_data -m paired --trim_primer --p1 <FORWARD_PRIMER> --p2 <REVERSE_PRIMER>
16s-preprocess -i ./raw_data -o ./cleaned_data -m single --skip_trim 
```

Tip: If primer detection rate is 0%, try swapping the sequences of `--p1` and `--p2`.

### 3. Chimera Removal

Remove potential chimeric sequences based on the reference database

``` bash
16s-rm-chimera -i ./cleaned_data -o ./no_chimera --db ref_database.db -k 25
```

### 4. Taxonomic Profiling

Generate genus-level abundance for each sample.

``` bash
16s-profile -i ./no_chimera -o ./abundance_out --db ref_database.db -k 25
```

### 5. Summary & Visualization

Merge all samples into a single matrix and generate composition plots.

``` bash
16s-summarize -i ./abundance_out -o ./final_report --top_n 20
```

## Expected Input Data

-   **Format**: Paired-end FASTQ files (`\_1.fastq.gz, \_2.fastq.gz`); Single-end FASTQ files (`*.fastq.gz`).

-   **Organization**: Each sample should ideally be in its own sub-directory or follow a consistent naming convention within the input folder.

-   **Current validation**: has been performed on oral microbiome 16S rRNA datasets using the file naming conventions described above. Additional FASTQ naming schemes may require minor adaptation.

## Key Outputs

-   `all_samples_genus_matrix.csv`: A unified matrix of genus abundances across all samples.
-   `composition.png`: A high-resolution stacked bar plot showing the top $N$ genera.
-   `*.log`: Detailed processing logs for reproducibility.

## Requirements

-   **Python**: 3.8+

-   **Core Dependencies**: Pandas, Numpy, Matplotlib, Seaborn, Scipy.

-   **Hardware**: Multi-threading is supported; 16GB+ RAM is recommended for large K-mer databases.

## Contact

#### Author: Xiaoqing Han

#### Email: [xhan723\@hotmail.com](mailto:xhan723@hotmail.com)
