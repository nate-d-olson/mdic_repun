# REPUN Pipeline

A Nextflow DSL2 pipeline for running [Repun](https://github.com/HKU-BAL/repun) (haplotype-aware variant representation unification) on multiple samples with support for BAM/CRAM files from local filesystem or S3.

## Features

- **Multi-sample processing** via samplesheet CSV
- **S3 and local file support** for BAM/CRAM inputs
- **Automatic CRAM to BAM conversion** with reference-based streaming
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

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV file |
| `--ref` | Local path to reference FASTA (with .fai index) |
| `--truth` | Local path to truth VCF file |

### Optional Arguments

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--cache_dir` | `${outdir}/s3_cache` | S3 file cache directory |
| `--force_download` | `false` | Force re-download of S3 files |

### Repun Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--somatic_mode` | `true` | Enable somatic variant mode |
| `--min_af` | `0.01` | Minimum allele frequency |
| `--max_af_somatic` | `0.25` | Max allele frequency for somatic unification |
| `--vaf_threshold` | `0.01` | VAF threshold for PASS variants |

### AWS Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aws_profile` | `mdic` | AWS profile name for S3 access |

## Samplesheet Format

Create a CSV file with the following columns:

```csv
sample_id,s3_bam_path,s3_bai_path,platform,output_subdir_name
HG002_ilmn,s3://bucket/HG002.bam,s3://bucket/HG002.bam.bai,ilmn,HG002_illumina_wgs
HG002_pacbio,s3://bucket/HG002.cram,s3://bucket/HG002.cram.crai,pacbio,HG002_pacbio_hifi
HG003_ont,/data/HG003.bam,/data/HG003.bam.bai,ont,HG003_nanopore_wgs
```

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier passed to Repun |
| `s3_bam_path` | S3 URI or local path to BAM/CRAM file |
| `s3_bai_path` | S3 URI or local path to index (.bai or .crai) |
| `platform` | Sequencing platform: `ilmn`, `pacbio`, or `ont` |
| `output_subdir_name` | Subdirectory name for outputs under `results/` |

**Notes:**
- CRAM files are automatically converted to BAM format
- S3 paths start with `s3://`, local paths are absolute file paths
- Index files (.crai) are converted to .bai for CRAM→BAM conversions

## Output Structure

Results are organized by sample under the output directory:

```
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

| Profile | CPUs | Memory | Use Case |
|---------|------|--------|----------|
| `standard` (default) | 4 | 32 GB | Default production runs |
| `local` | 24 | 386 GB | High-resource local machine |
| `high_mem` | 22 | 128 GB | Memory-intensive samples |
| `test` | 2 | 6 GB | Quick testing with minimal resources |
| `debug` | - | - | Debug mode with cleanup disabled |
| `conda` | - | - | Use Conda instead of Docker |

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

## AWS Configuration

Ensure AWS credentials are available via one of:

1. **Environment variables:**
   ```bash
   export AWS_ACCESS_KEY_ID=your_access_key
   export AWS_SECRET_ACCESS_KEY=your_secret_key
   ```

2. **AWS credentials file:** `~/.aws/credentials`
   ```ini
   [mdic]
   aws_access_key_id = your_access_key
   aws_secret_access_key = your_secret_key
   ```

3. **IAM instance role** (if running on EC2)

The pipeline uses the `mdic` profile by default. Change this via `--aws_profile`.

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

## Platform-Specific Notes

### Illumina (`ilmn`)
- Short-read sequencing data
- Typically WGS or WES BAM files
- Works with both BAM and CRAM formats

### PacBio HiFi (`pacbio`)
- Long-read HiFi sequencing
- High accuracy reads
- CRAM format commonly used

### Oxford Nanopore (`ont`)
- Long-read sequencing
- Handles both R9 and R10 chemistry
- BAM/CRAM format supported

## Citation

If you use this pipeline, please cite:

**Repun:**
> Zheng Z, Li S, Su J, et al. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nature Computational Science* 2, 797–803 (2022). https://doi.org/10.1038/s43588-022-00387-x

**Nextflow:**
> Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35, 316–319 (2017). https://doi.org/10.1038/nbt.3820

## License

This pipeline wrapper is provided as-is. See the [Repun repository](https://github.com/HKU-BAL/repun) for Repun-specific licensing.
