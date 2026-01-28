# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a **Nextflow DSL2 pipeline wrapper** for running [Repun](https://github.com/HKU-BAL/repun) on multiple samples with S3/CRAM input support. The repository contains both the Nextflow orchestration layer and a clone of the Repun tool itself.

**Repun** is a haplotype-aware variant representation unification algorithm that harmonizes variant representations between truth VCF files and BAM/CRAM alignments across different sequencing platforms (ONT, PacBio HIFI, Illumina).

## Directory Structure

```
/wrk/mdic_repun/
├── main.nf                    # Nextflow DSL2 pipeline entry point
├── nextflow.config            # Nextflow configuration and profiles
├── README.md                  # User documentation
├── CLAUDE.md                  # AI assistant guidance (this file)
├── samplesheet*.csv          # Sample input configurations
├── Repun/                    # Repun tool submodule (see Repun/CLAUDE.md)
├── data/                     # Local data directory
├── results/                  # Pipeline output directory
└── work/                     # Nextflow work directory (intermediate files)
```

## Running the Pipeline

The pipeline orchestrates Repun execution across multiple samples defined in a samplesheet:

```bash
# Basic execution
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

# Disable somatic mode
nextflow run main.nf \
    --input samplesheet.csv \
    --ref /path/to/reference.fasta \
    --truth /path/to/truth.vcf \
    --somatic_mode false
```

## Nextflow Configuration Profiles

Defined in [nextflow.config](nextflow.config):

| Profile | CPUs | Memory | Use Case |
|---------|------|--------|----------|
| **standard** (default) | 4 | 32 GB | Default production runs |
| **local** | 24 | 386 GB | High-resource local machine |
| **high_mem** | 22 | 128 GB | Memory-intensive samples |
| **test** | 2 | 6 GB | Quick testing with minimal resources |
| **debug** | - | - | Debug mode with cleanup disabled |
| **conda** | - | - | Use Conda instead of Docker |

## Samplesheet Format

The pipeline expects a CSV with the following columns:

```csv
sample_id,s3_bam_path,s3_bai_path,platform,output_subdir_name
```

- `sample_id`: Sample identifier passed to Repun's `--sample` parameter
- `s3_bam_path`: S3 URI to BAM or CRAM file (supports local paths too)
- `s3_bai_path`: S3 URI to BAM/CRAM index (.bai, .crai)
- `platform`: Sequencing platform - must be `ilmn`, `pacbio`, or `ont`
- `output_subdir_name`: Subdirectory name under `results/` for this sample's output

## Pipeline Parameters

### Required Parameters
- `--input`: Path to samplesheet CSV
- `--ref`: Local path to reference FASTA (must have .fai index)
- `--truth`: Local path to truth VCF file

### Optional Parameters
- `--outdir`: Output directory (default: `results`)
- `--cache_dir`: S3 file cache directory (default: `${outdir}/s3_cache`)
- `--force_download`: Force re-download of S3 files (default: `false`)

### Repun Configuration Parameters
- `--somatic_mode`: Enable somatic variant mode (default: `true`)
- `--min_af`: Minimum allele frequency (default: `0.01`)
- `--max_af_somatic`: Max allele frequency for somatic unification (default: `0.25`)
- `--vaf_threshold`: VAF threshold for PASS variants (default: `0.01`)

### AWS Configuration Parameters
- `--aws_profile`: AWS profile name for S3 access (default: `mdic`)

## AWS Configuration

The pipeline uses the `mdic` AWS profile by default (configured in [nextflow.config:90-100](nextflow.config)). Ensure credentials are available via:

- AWS credentials file: `~/.aws/credentials` with `[mdic]` profile
- Environment variables: `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`
- IAM instance role (if running on EC2)

## Reference Files

The pipeline requires:
- Reference FASTA file (specified via `--ref`)
- Reference FASTA index (.fai file, must exist alongside the FASTA)
- Truth VCF file (specified via `--truth`)

For CRAM files, ensure reference genome is accessible. The pipeline sets appropriate environment variables for CRAM reference cache automatically.

## Repun Tool Details

For detailed information about the Repun tool itself, see `Repun/CLAUDE.md`. Key points:

### Running Repun Directly

```bash
# Inside Docker container
docker run -it \
  -v /data:/data \
  hkubal/repun:latest \
  /opt/bin/repun \
    --bam_fn /data/sample.bam \
    --ref_fn /data/ref.fa \
    --truth_vcf_fn /data/truth.vcf \
    --threads 24 \
    --platform ilmn \
    --output_dir /data/output

# Somatic mode
docker run -it \
  -v /data:/data \
  hkubal/repun:latest \
  /opt/bin/repun \
    --bam_fn /data/sample.bam \
    --ref_fn /data/ref.fa \
    --truth_vcf_fn /data/truth.vcf \
    --threads 24 \
    --platform ilmn \
    --somatic_mode \
    --max_af_for_somatic_unification 0.01 \
    --vaf_threshold_for_pass 0.01 \
    --output_dir /data/output
```

### Platforms
- `ilmn`: Illumina short-read sequencing
- `hifi` or `pacbio`: PacBio HIFI long-read sequencing
- `ont`: Oxford Nanopore Technology long-read sequencing

## Output Structure

```
results/
├── <output_subdir_name_1>/
│   └── repun_output/
│       ├── unified.vcf.gz           # Germline mode output
│       ├── output_truth.vcf.gz      # Somatic mode (truth coordinates)
│       └── output_candidate.vcf.gz  # Somatic mode (candidate coordinates)
├── <output_subdir_name_2>/
│   └── repun_output/
│       └── ...
├── pipeline_info/
│   ├── execution_timeline.html      # Execution timeline
│   ├── execution_report.html        # Resource usage report
│   ├── execution_trace.txt          # Detailed trace file
│   └── pipeline_dag.svg             # Pipeline DAG visualization
└── s3_cache/                        # Cached S3 downloads (if using S3)
    └── s3_downloads/
```

## Important Implementation Notes

### Container Configuration
The pipeline is currently configured to use a Conda environment for development:
- Current: Uses hardcoded Conda environment path in [main.nf:104](main.nf)
- Production: Will use `hkubal/repun:latest` Docker container (currently commented out)
- Note: Conda configuration is being used while optimizing CRAM file handling performance

### S3 File Caching
The pipeline intelligently caches S3 downloads to avoid redundant transfers:
- Downloaded files are stored in `${params.cache_dir}/s3_downloads/`
- Reused across pipeline runs unless `--force_download` is specified
- Helps reduce data transfer costs and improves performance

### Nextflow Work Directory
Nextflow stores intermediate files and logs in `work/`.
- Use `-resume` to reuse cached results from previous runs
- Clean with `nextflow clean -f` to free disk space
- Check `.nextflow.log` for detailed execution logs

## Common Development Workflows

### Customizing Repun Parameters
Repun parameters can be customized via command-line options:

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

### Using Different Execution Profiles

```bash
# High-resource local machine
nextflow run main.nf -profile local --input samplesheet.csv --ref ref.fa --truth truth.vcf

# Memory-intensive samples
nextflow run main.nf -profile high_mem --input samplesheet.csv --ref ref.fa --truth truth.vcf

# Testing with minimal resources
nextflow run main.nf -profile test --input samplesheet.csv --ref ref.fa --truth truth.vcf
```

### Debugging Failed Runs
```bash
# Check Nextflow logs
cat .nextflow.log

# Inspect work directory for specific task
ls -la work/<task-hash>/

# Resume from checkpoint
nextflow run main.nf -resume --input samplesheet.csv --ref ref.fa --truth truth.vcf
```

## Dependencies

- Nextflow >= 21.04
- Docker (recommended) or Conda/Mamba
- AWS CLI configured (for S3 access)
- Sufficient disk space in work directory for intermediate files

## Pipeline Reports

After execution, comprehensive reports are generated in `results/pipeline_info/`:
- **execution_report.html**: Resource usage and task statistics
- **execution_timeline.html**: Visual timeline of task execution
- **execution_trace.txt**: Detailed trace of all tasks
- **pipeline_dag.svg**: Pipeline dependency graph