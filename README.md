# REPUN Pipeline

A simple Nextflow DSL2 pipeline for running [repun](https://github.com/HKU-BAL/repun) on multiple samples with S3 input support.

## Requirements

- Nextflow >= 21.04
- Docker
- AWS CLI configured (for S3 access)

## Quick Start

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref /path/to/reference.fasta \
    --truth /path/to/truth.vcf
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--input` | Yes | - | Path to samplesheet CSV |
| `--ref` | Yes | - | Local path to reference FASTA (with .fai index) |
| `--truth` | Yes | - | Local path to truth VCF |
| `--outdir` | No | `results` | Output directory |

## Samplesheet Format

Create a CSV file with the following columns:

```csv
sample_id,s3_bam_path,s3_bai_path,platform,output_subdir_name
broad_pacbio_40x,s3://bucket/sample.bam,s3://bucket/sample.bam.bai,pacbio,run_name
```

| Column | Description |
|--------|-------------|
| `sample_id` | Sample name passed to repun `--sample` |
| `s3_bam_path` | S3 URI to BAM or CRAM file |
| `s3_bai_path` | S3 URI to BAM/CRAM index |
| `platform` | Sequencing platform: `ilmn` or `pacbio` |
| `output_subdir_name` | Subdirectory name for outputs under `results/` |

## Output

Results are organized by sample under the output directory:

```
results/
├── broad_pacbio_40x_wgs/
│   └── <repun output files>
├── broad_ilmn_120x_wgs/
│   └── <repun output files>
└── ffpe_rna_illumina/
    └── <repun output files>
```

## Profiles

| Profile | CPUs | Memory | Use Case |
|---------|------|--------|----------|
| `standard` | 22 | 88 GB | Default production runs |
| `local` | 8 | 32 GB | Local testing |
| `high_mem` | 22 | 128 GB | Memory-intensive samples |

```bash
nextflow run main.nf -profile local --input samplesheet.csv --ref ref.fa --truth truth.vcf
```

## AWS Configuration

Ensure AWS credentials are available via one of:

- Environment variables: `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`
- AWS credentials file: `~/.aws/credentials`
- IAM instance role (if running on EC2)

## Resource Configuration

Default resources (22 CPUs, 88 GB memory) can be overridden in `nextflow.config` or via command line:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --ref ref.fa \
    --truth truth.vcf \
    -process.cpus 16 \
    -process.memory '64 GB'
```

## Resume Failed Runs

Nextflow caches completed tasks. Resume a failed run with:

```bash
nextflow run main.nf -resume --input samplesheet.csv --ref ref.fa --truth truth.vcf
```
