# REPUN Pipeline

A Nextflow DSL2 pipeline for running [Repun](https://github.com/HKU-BAL/repun) (haplotype-aware variant representation unification) on multiple samples with support for BAM/CRAM files from local filesystem or S3.

## Features

- **Multi-sample processing** via `samplesheet.csv`
- **S3 and local file support** for BAM/CRAM inputs
- **Intelligent file caching** to avoid redundant downloads
- **Configurable somatic mode** with tunable parameters
- **Multiple execution profiles** for different compute environments
- **Comprehensive reporting** with timeline, trace, and DAG visualization

## Requirements

- Nextflow >= 21.04
- Docker (recommended) or Conda/Mamba
- AWS CLI configured (for S3 access)

## Quick Start

```bash
# Basic run
nextflow run main.nf \
    --input samplesheet.csv \
    --ref /path/to/reference.fasta \
    --truth /path/to/truth.vcf

# With specific profile
nextflow run main.nf -profile local \
    --input samplesheet.csv \
    --ref /path/to/reference.fasta \
    --truth /path/to/truth.vcf

# Resume failed run
nextflow run main.nf -resume \
    --input samplesheet.csv \
    --ref /path/to/reference.fasta \
    --truth /path/to/truth.vcf
```

## Parameters

### Required Arguments

| Parameter | Description                                     |
| --------- | ----------------------------------------------- |
| `--input` | Path to samplesheet CSV file                    |
| `--ref`   | Local path to reference FASTA (with .fai index) |
| `--truth` | Local path to truth VCF file                    |

### Optional Arguments

| Parameter          | Default              | Description                   |
| ------------------ | -------------------- | ----------------------------- |
| `--outdir`         | `results`            | Output directory              |
| `--cache_dir`      | `${outdir}/s3_cache` | S3 file cache directory       |
| `--force_download` | `false`              | Force re-download of S3 files |

### Repun Configuration

| Parameter          | Default | Description                                  |
| ------------------ | ------- | -------------------------------------------- |
| `--somatic_mode`   | `true`  | Enable somatic variant mode                  |
| `--min_af`         | `0.01`  | Minimum allele frequency                     |
| `--max_af_somatic` | `0.25`  | Max allele frequency for somatic unification |
| `--vaf_threshold`  | `0.01`  | VAF threshold for PASS variants              |

### AWS Configuration

| Parameter       | Default | Description                    |
| --------------- | ------- | ------------------------------ |
| `--aws_profile` | `mdic`  | AWS profile name for S3 access |

Sample sheets included in this repo are specifically for the MDIC SRS project.
These files are not public and provided for transparency, not intended for use 
outside of the NIST-MDIC SRS bioinformatics team.

## Samplesheet Format

Create a CSV file with the following columns:

```csv
sample_id,s3_bam_path,s3_bai_path,platform,output_subdir_name
HG002_ilmn,s3://bucket/HG002.bam,s3://bucket/HG002.bam.bai,ilmn,HG002_illumina_wgs
HG002_pacbio,s3://bucket/HG002.cram,s3://bucket/HG002.cram.crai,pacbio,HG002_pacbio_hifi
HG003_ont,/data/HG003.bam,/data/HG003.bam.bai,ont,HG003_nanopore_wgs
```

| Column               | Description                                     |
| -------------------- | ----------------------------------------------- |
| `sample_id`          | Sample identifier passed to Repun               |
| `s3_bam_path`        | S3 URI or local path to BAM/CRAM file           |
| `s3_bai_path`        | S3 URI or local path to index (.bai or .crai)   |
| `platform`           | Sequencing platform: `ilmn`, `pacbio`, or `ont` |
| `output_subdir_name` | Subdirectory name for outputs under `results/`  |

## Output Structure

Results are organized by sample under the output directory:

```txt
results/
├── HG002_illumina_wgs/
│   └── repun_output/
│       ├── unified.vcf.gz              # Germline mode output
│       ├── output_truth.vcf.gz         # Somatic mode (truth coordinates)
│       └── output_candidate.vcf.gz     # Somatic mode (candidate coordinates)
├── HG002_pacbio_hifi/
│   └── repun_output/
│       └── ...
├── pipeline_info/
│   ├── execution_timeline.html         # Execution timeline
│   ├── execution_report.html           # Resource usage report
│   ├── execution_trace.txt             # Detailed trace file
│   └── pipeline_dag.svg                # Pipeline DAG visualization
└── s3_cache/                           # Cached S3 downloads
    └── s3_downloads/
        ├── sample1.bam
        └── sample1.bam.bai
```

## Execution Profiles

The pipeline includes several pre-configured profiles for different compute environments:

| Profile              | CPUs | Memory | Use Case                             |
| -------------------- | ---- | ------ | ------------------------------------ |
| `standard` (default) | 4    | 32 GB  | Default production runs              |
| `local`              | 24   | 386 GB | High-resource local machine          |
| `high_mem`           | 22   | 128 GB | Memory-intensive samples             |
| `test`               | 2    | 6 GB   | Quick testing with minimal resources |
| `debug`              | -    | -      | Debug mode with cleanup disabled     |
| `conda`              | -    | -      | Use Conda instead of Docker          |

### Using Profiles

```bash
# Local high-resource execution
nextflow run main.nf -profile local \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf

# Memory-intensive samples
nextflow run main.nf -profile high_mem \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf

# Debug mode
nextflow run main.nf -profile debug \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf
```

## Advanced Usage

### Disable Somatic Mode

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf \
    --somatic_mode false
```

### Custom Somatic Parameters

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf \
    --somatic_mode true \
    --min_af 0.05 \
    --max_af_somatic 0.3 \
    --vaf_threshold 0.05
```

### Override Resource Limits

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf \
    --max_cpus 32 \
    --max_memory '512.GB' \
    --max_time '48.h'
```

### Use Custom AWS Profile

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf \
    --aws_profile my_custom_profile
```

## Troubleshooting

### Resume Failed Runs

Nextflow caches completed tasks. Resume a failed run with:

```bash
nextflow run main.nf -resume \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf
```

### View Help Message

```bash
nextflow run main.nf --help
```

### Check Pipeline Reports

After execution, check the reports in `results/pipeline_info/`:

- `execution_report.html` - Resource usage and task statistics
- `execution_timeline.html` - Visual timeline of task execution
- `execution_trace.txt` - Detailed trace of all tasks
- `pipeline_dag.svg` - Pipeline dependency graph

### Clean Work Directory

Remove intermediate files to free disk space:

```bash
nextflow clean -f
```

### Debug Mode

Run in debug mode to preserve work directories:

```bash
nextflow run main.nf -profile debug \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf
```

## Citation

If you use this pipeline, please cite:

**Repun:**

> Zheng Z, Li S, Su J, et al. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. _Nature Computational Science_ 2, 797–803 (2022). <https://doi.org/10.1038/s43588-022-00387-x>

**Nextflow:**

> Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. _Nat Biotechnol_ 35, 316–319 (2017). <https://doi.org/10.1038/nbt.3820>

## AI Use Disclosure

This pipeline was developed with assistance from Anthropic's Claude Code (Claude Sonnet 4.5) for pipeline development, code generation, optimization, and documentation. All AI-generated code has been manually tested and validated by NIST researchers. No restricted access data was provided to external AI services during development.

## License

This data/work was created by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States. This data/work may be subject to foreign copyright.

The data/work is provided by NIST as a public service and is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NIST does not warrant or make any representations regarding the use of the data or the results thereof, including but not limited to the correctness, accuracy, reliability or usefulness of the data. NIST SHALL NOT BE LIABLE AND YOU HEREBY RELEASE NIST FROM LIABILITY FOR ANY INDIRECT, CONSEQUENTIAL, SPECIAL, OR INCIDENTAL DAMAGES (INCLUDING DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION, LOSS OF BUSINESS INFORMATION, AND THE LIKE), WHETHER ARISING IN TORT, CONTRACT, OR OTHERWISE, ARISING FROM OR RELATING TO THE DATA (OR THE USE OF OR INABILITY TO USE THIS DATA), EVEN IF NIST HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

To the extent that NIST may hold copyright in countries other than the United States, you are hereby granted the non-exclusive irrevocable and unconditional right to print, publish, prepare derivative works and distribute the NIST data, in any medium, or authorize others to do so on your behalf, on a royalty-free basis throughout the world.

You may improve, modify, and create derivative works of the data or any portion of the data, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the data and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the data: Data citation recommendations are provided at <https://www.nist.gov/open/license>.

Permission to use this data is contingent upon your acceptance of the terms of this agreement and upon your providing appropriate acknowledgments of NIST's creation of the data/work.
