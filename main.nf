#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    REPUN Pipeline
========================================================================================
    A Nextflow DSL2 pipeline for running Repun (haplotype-aware variant representation
    unification) on multiple samples with support for BAM/CRAM files from local or S3.

    GitHub: https://github.com/HKU-BAL/repun
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

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
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
// Process: DOWNLOAD_AND_CONVERT_ALIGNMENT
// =============================================================================

process DOWNLOAD_AND_CONVERT_ALIGNMENT {
    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::samtools=1.17'

    publishDir "${params.cache_dir}/s3_downloads",
        mode: 'copy',
        overwrite: params.force_download,
        enabled: params.cache_dir ? true : false

    input:
    tuple val(sample_id), val(file_path), val(platform), val(output_name)
    path ref
    path ref_fai

    output:
    tuple val(sample_id), path("${output_file_name}"), val(platform), val(output_name)

    script:
    is_s3 = file_path.startsWith('s3://')
    is_cram = file_path.endsWith('.cram')
    base_name = file_path.tokenize('/').last()

    // Output is always .bam (converted from .cram if needed)
    output_file_name = is_cram ? base_name.replaceAll('\\.cram$', '.bam') : base_name

    if (is_cram) {
        """
        set -e
        echo "[CONVERT] Streaming and converting CRAM to BAM: ${file_path}"
        echo "Sample: ${sample_id}"
        echo "Platform: ${platform}"

        samtools version

        # Stream CRAM from S3/local and convert to BAM
        export AWS_PROFILE=${params.aws_profile}

        samtools view \\
            -b \\
            -T ${ref} \\
            -@ ${task.cpus} \\
            -o ${output_file_name} \\
            ${file_path}

        if [[ ! -f "${output_file_name}" ]]; then
            echo "ERROR: Failed to convert CRAM ${file_path}"
            exit 1
        fi

        echo "[CONVERT] Successfully converted to BAM"
        ls -lh "${output_file_name}"
        """
    } else if (is_s3) {
        """
        set -e
        echo "[DOWNLOAD] Downloading BAM from S3: ${file_path}"
        export AWS_PROFILE=${params.aws_profile}

        aws s3 cp ${file_path} ${output_file_name}

        if [[ ! -f "${output_file_name}" ]]; then
            echo "ERROR: Failed to download ${file_path}"
            exit 1
        fi

        echo "[DOWNLOAD] Successfully downloaded BAM"
        ls -lh "${output_file_name}"
        """
    } else {
        """
        echo "[LOCAL] Using local file: ${file_path}"
        ln -s ${file_path} ${output_file_name}
        ls -lh "${output_file_name}"
        """
    }
}

// =============================================================================
// Process: HANDLE_INDEX
// =============================================================================

process HANDLE_INDEX {
    tag "${sample_id}"
    label 'process_low'

    conda 'bioconda::samtools=1.17'

    publishDir "${params.cache_dir}/s3_downloads",
        mode: 'copy',
        overwrite: params.force_download,
        enabled: params.cache_dir ? true : false

    input:
    tuple val(sample_id), val(index_path), val(platform), val(output_name)
    tuple val(sample_id), path(bam_file), val(bam_platform), val(bam_output_name)

    output:
    tuple val(sample_id), path("${output_file_name}"), val(platform), val(output_name)

    script:
    is_s3 = index_path.startsWith('s3://')
    is_crai = index_path.endsWith('.crai')
    base_name = index_path.tokenize('/').last()

    // Output is always .bai (generated from BAM if source was .crai)
    output_file_name = is_crai ? base_name.replaceAll('\\.crai$', '.bai') : base_name

    if (is_crai) {
        """
        set -e
        echo "[INDEX] Generating BAI index from converted BAM: ${bam_file}"
        echo "Sample: ${sample_id}"

        samtools index -@ ${task.cpus} ${bam_file} ${output_file_name}

        if [[ ! -f "${output_file_name}" ]]; then
            echo "ERROR: Failed to generate index for ${bam_file}"
            exit 1
        fi

        echo "[INDEX] Successfully created BAI index"
        ls -lh "${output_file_name}"
        """
    } else if (is_s3) {
        """
        set -e
        echo "[DOWNLOAD] Downloading index from S3: ${index_path}"
        export AWS_PROFILE=${params.aws_profile}

        aws s3 cp ${index_path} ${output_file_name}

        if [[ ! -f "${output_file_name}" ]]; then
            echo "ERROR: Failed to download ${index_path}"
            exit 1
        fi

        echo "[DOWNLOAD] Successfully downloaded index"
        ls -lh "${output_file_name}"
        """
    } else {
        """
        echo "[LOCAL] Using local index: ${index_path}"
        ln -s ${index_path} ${output_file_name}
        ls -lh "${output_file_name}"
        """
    }
}

// =============================================================================
// Workflow
// =============================================================================

workflow {

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
    Channel
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

    // Create separate channels for alignment and index files
    alignment_files_ch = samples_raw_ch.map { sample_id, bam_path, bai_path, platform, output_name ->
        tuple(sample_id, bam_path, platform, output_name)
    }

    index_files_ch = samples_raw_ch.map { sample_id, bam_path, bai_path, platform, output_name ->
        tuple(sample_id, bai_path, platform, output_name)
    }

    // Download/convert alignment files (BAM or CRAM → BAM)
    // DOWNLOAD_AND_CONVERT_ALIGNMENT(
    //     alignment_files_ch,
    //     ref_file,
    //     ref_fai_file
    // )

    // Handle index files (download .bai or generate from converted BAM if source was .crai)
    // Join with alignment files so we have the BAM available for indexing
    // index_files_ch
    //     .join(DOWNLOAD_AND_CONVERT_ALIGNMENT.out, by: 0)  // Join by sample_id
    //     .multiMap { sample_id, idx_path, idx_platform, idx_output, bam_file, bam_platform, bam_output ->
    //         index_input: tuple(sample_id, idx_path, idx_platform, idx_output)
    //         bam_input: tuple(sample_id, bam_file, bam_platform, bam_output)
    //     }
    //     .set { index_inputs }

    // HANDLE_INDEX(
    //     index_inputs.index_input,
    //     index_inputs.bam_input
    // )

    // // Combine alignment and index files for RUN_REPUN
    // DOWNLOAD_AND_CONVERT_ALIGNMENT.out
    //     .join(HANDLE_INDEX.out, by: 0)  // Join by sample_id
    //     .map { sample_id, bam_file, bam_platform, bam_output, bai_file, bai_platform, bai_output ->
    //         tuple(sample_id, bam_file, bai_file, bam_platform, bam_output)
    //     }
    //     .set { samples_prepared_ch }

    // Run repun on each sample
    RUN_REPUN(
        samples_raw_ch,
        ref_file,
        ref_fai_file,
        truth_file
    )
}
