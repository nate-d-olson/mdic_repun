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
        --local_data_dir <path>  Local directory to check for cached files (default: data)
                                 Set to 'false' to disable local file detection

    Repun Options:
        --somatic_mode       Enable somatic mode (default: true)
        --min_af <float>     Minimum allele frequency (default: 0.01)
        --max_af_somatic <float>  Max AF for somatic unification (default: 0.01)
        --vaf_threshold <float>   VAF threshold for PASS (default: 0.01)
        --chunk_size <int>		Size of chunks to split contigs into for parallelization

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

        # Force S3 download (ignore local files)
        nextflow run main.nf --input samples.csv --ref ref.fa --truth truth.vcf --local_data_dir false

    For more information, visit: https://github.com/HKU-BAL/repun
    """.stripIndent()
}

// =============================================================================
// Helper Functions
// =============================================================================

/**
 * Resolve file path: check if local version exists, otherwise use remote path.
 * Searches for file by basename in the local data directory.
 *
 * @param remote_path The S3 or remote file path from samplesheet
 * @param local_dir   The local directory to search for cached files
 * @return The local path if file exists locally, otherwise the original remote path
 */
def resolveFilePath(remote_path, local_dir) {
    if (!local_dir || local_dir == 'false' || local_dir == false) {
        return remote_path
    }

    def basename = file(remote_path).getName()
    def local_path = "${local_dir}/${basename}"
    def local_file = file(local_path)

    if (local_file.exists()) {
        log.info "  Using local file: ${local_path}"
        return local_path
    } else {
        log.info "  Using remote file: ${remote_path}"
        return remote_path
    }
}

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
    path truth_tbi

    output:
    path "repun_*/*", emit: results
    path ".command.log", emit: log, optional: true

    script:
    def somatic_args = params.somatic_mode ? """\\
        --somatic_mode \\
        --min_af ${params.min_af} \\
        --max_af_for_somatic_unification ${params.max_af_somatic} \\
        --vaf_threshold_for_pass ${params.vaf_threshold} \\
        --chunk_size ${params.chunk_size}""" : ""

    def output_subdir = params.somatic_mode
        ? "repun_sm_af${params.min_af}_maxaf${params.max_af_somatic}_vaf${params.vaf_threshold}"
        : "repun_germline"

    """
    mkdir -p ${output_subdir}

    echo "Starting Repun for sample: ${sample_id}"
    echo "Platform: ${platform}"
    echo "Somatic mode: ${params.somatic_mode}"
    echo "Output directory: ${output_subdir}"

    Repun/repun \\
        --bam_fn ${bam} \\
        --ref_fn ${ref} \\
        --truth_vcf_fn ${truth} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir ${output_subdir} \\
        --sample_name ${sample_id} \\
        ${somatic_args}

    echo "Repun completed for sample: ${sample_id}"
    ls -lh ${output_subdir}/
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
    Local data dir    : ${params.local_data_dir ?: 'disabled'}
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
    truth_tbi_file = file("${params.truth}.tbi", checkIfExists: true)


    // Resolve local data directory path
    def local_data_dir = params.local_data_dir
    if (local_data_dir && local_data_dir != 'false' && local_data_dir != false) {
        def local_dir_file = file(local_data_dir)
        if (!local_dir_file.exists()) {
            log.warn "Local data directory does not exist: ${local_data_dir} - will use remote paths"
            local_data_dir = false
        } else {
            log.info "Checking for local files in: ${local_data_dir}"
        }
    }

    // Parse samplesheet and check for local file versions
    channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id   = row.sample_id.trim()
            def bam_path    = resolveFilePath(row.s3_bam_path.trim(), local_data_dir)
            def bai_path    = resolveFilePath(row.s3_bai_path.trim(), local_data_dir)
            def platform    = row.platform.trim()
            def output_name = row.output_subdir_name.trim()

            tuple(sample_id, bam_path, bai_path, platform, output_name)
        }
        .set { samples_raw_ch }


    // Run repun on each sample
    RUN_REPUN(
        samples_raw_ch,
        ref_file,
        ref_fai_file,
        truth_file,
        truth_tbi_file
    )
}
