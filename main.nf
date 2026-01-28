#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    REPUN Pipeline
========================================================================================
    A Nextflow DSL2 pipeline for running Repun (haplotype-aware variant representation
    unification) on multiple samples with support for BAM/CRAM files from local or S3.

    GitHub: https://github.com/HKU-BAL/repun

    AI Use Disclosure:
    This pipeline utilized Anthropic's Claude Code (Claude Sonnet 4.5) for pipeline
    development, code generation, and documentation. All code has been reviewed and
    validated by NIST researchers.
----------------------------------------------------------------------------------------
*/

// =============================================================================
// Help Message
// =============================================================================

def helpMessage() {
    log.info"""
    Usage:
        nextflow run main.nf [options]

    Required Arguments:
        --input <path>       Path to samplesheet CSV file
        --ref <path>         Local path to reference FASTA file (with .fai index)
        --truth <path>       Local path to truth VCF file

    Optional Arguments:
        --outdir <path>      Output directory (default: results)
        --cache_dir <path>   S3 file cache directory (default: \${outdir}/s3_cache)
        --force_download     Force re-download of S3 files (default: false)

    Repun Options:
        --somatic_mode       Enable somatic mode (default: true)
        --min_af <float>     Minimum allele frequency (default: 0.01)
        --max_af_somatic <float>  Max AF for somatic unification (default: 0.25)
        --vaf_threshold <float>   VAF threshold for PASS (default: 0.01)

    AWS Options:
        --aws_profile <str>  AWS profile name (default: mdic)

    Profiles:
        -profile standard    Default - 4 CPUs, 32 GB memory
        -profile local       High resource - 24 CPUs, 386 GB memory
        -profile high_mem    Memory intensive - 22 CPUs, 128 GB memory

    Examples:
        # Basic run
        nextflow run main.nf --input samples.csv --ref ref.fa --truth truth.vcf

        # With specific profile
        nextflow run main.nf -profile local --input samples.csv --ref ref.fa --truth truth.vcf

        # Disable somatic mode
        nextflow run main.nf --input samples.csv --ref ref.fa --truth truth.vcf --somatic_mode false

    For more information, visit: https://github.com/HKU-BAL/repun
    """.stripIndent()
}

// =============================================================================
// Parameters
// =============================================================================

// Required inputs
params.input  = null   // Path to samplesheet CSV
params.ref    = null   // Local path to reference FASTA
params.truth  = null   // Local path to truth VCF

// Output configuration
params.outdir = 'results'

// S3 file caching configuration
params.cache_dir = "${params.outdir}/s3_cache"
params.force_download = false  // Set to true to force re-download

// Repun somatic mode parameters
params.somatic_mode = true
params.min_af = 0.01
params.max_af_somatic = 0.25
params.vaf_threshold = 0.01

// AWS configuration
params.aws_profile = 'mdic'

// Help parameter
params.help = false

// =============================================================================
// Process: RUN_REPUN
// =============================================================================

process RUN_REPUN {
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/${output_name}",
        mode: 'copy',
        saveAs: { filename -> filename.equals('.command.log') ? null : filename }

    // container 'hkubal/repun:latest'
    conda '/home/nolson/miniforge/envs/run'

    input:
    tuple val(sample_id), path(bam), path(bai), val(platform), val(output_name)
    path ref
    path ref_fai
    path truth

    output:
    path "repun_output/*", emit: results
    path ".command.log", emit: log, optional: true

    script:
    def somatic_args = params.somatic_mode ? """\\
        --somatic_mode \\
        --min_af ${params.min_af} \\
        --max_af_for_somatic_unification ${params.max_af_somatic} \\
        --vaf_threshold_for_pass ${params.vaf_threshold}""" : ""

    """
    mkdir -p repun_output

    echo "Starting Repun for sample: ${sample_id}"
    echo "Platform: ${platform}"
    echo "Somatic mode: ${params.somatic_mode}"

	vmtouch -dl ${bai} ${ref} ${ref_fai}
	
	export REF_PATH=/resources/ref_cache/%2s/%2s/%s
	export REF_CACHE=/resources/ref_cache/%2s/%2s/%s
	
    python /wrk/mdic_repun/Repun/repun \\
        --bam_fn ${bam} \\
        --ref_fn ${ref} \\
        --truth_vcf_fn ${truth} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir repun_output ${somatic_args}

    echo "Repun completed for sample: ${sample_id}"
    ls -lh repun_output/
    """
}


// =============================================================================
// Workflow
// =============================================================================

workflow {

    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Print pipeline header
    log.info """\
    ============================================
    R E P U N   P I P E L I N E
    ============================================
    Input samplesheet : ${params.input}
    Reference FASTA   : ${params.ref}
    Truth VCF         : ${params.truth}
    Output directory  : ${params.outdir}
    Cache directory   : ${params.cache_dir}
    Somatic mode      : ${params.somatic_mode}
    AWS profile       : ${params.aws_profile}
    ============================================
    """.stripIndent()

    // Parameter Validation
    if (!params.input) {
        error "ERROR: Please provide a samplesheet via --input"
    }
    if (!params.ref) {
        error "ERROR: Please provide reference FASTA via --ref"
    }
    if (!params.truth) {
        error "ERROR: Please provide truth VCF via --truth"
    }

    // Load reference files early (needed for CRAM conversion)
    ref_file     = file(params.ref, checkIfExists: true)
    ref_fai_file = file("${params.ref}.fai", checkIfExists: true)
    truth_file   = file(params.truth, checkIfExists: true)


    // Parse samplesheet and separate alignment files from index files
    channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id   = row.sample_id
            def bam_path    = row.s3_bam_path
            def bai_path    = row.s3_bai_path
            def platform    = row.platform
            def output_name = row.output_subdir_name

            tuple(sample_id, bam_path, bai_path, platform, output_name)
        }
        .set { samples_raw_ch }


    // Run repun on each sample
    RUN_REPUN(
        samples_raw_ch,
        ref_file,
        ref_fai_file,
        truth_file
    )
}
