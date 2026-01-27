#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * REPUN Pipeline
 * A simple DSL2 pipeline to run repun on multiple samples from S3
 */

// =============================================================================
// Parameters
// =============================================================================

params.input  = null   // Path to samplesheet CSV
params.ref    = null   // Local path to reference FASTA
params.truth  = null   // Local path to truth VCF
params.outdir = 'results'

// S3 file caching configuration
params.cache_dir = "${params.outdir}/s3_cache"
params.force_download = false  // Set to true to force re-download

// =============================================================================
// Process: RUN_REPUN
// =============================================================================

process RUN_REPUN {
    tag "${sample_id}"
    publishDir "${params.outdir}/${output_name}", mode: 'copy'

    container 'hkubal/repun:latest'

    input:
    tuple val(sample_id), path(bam), path(bai), val(platform), val(output_name)
    path ref
    path ref_fai
    path truth

    output:
    path "repun_output/*", emit: results

    script:
    """
    mkdir -p repun_output

    repun \\
        --bam_fn ${bam} \\
        --ref_fn ${ref} \\
        --truth_vcf_fn ${truth} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir repun_output \\
        --somatic_mode \\
        --min_af 0.01 \\
        --max_af_for_somatic_unification 0.25 \\
        --vaf_threshold_for_pass 0.01
    """
}

// =============================================================================
// Process: DOWNLOAD_AND_CONVERT_ALIGNMENT
// =============================================================================

process DOWNLOAD_AND_CONVERT_ALIGNMENT {
    tag "${sample_id}"

    conda 'samtools'

    // Publish to cache
    publishDir "${params.cache_dir}/s3_downloads",
        mode: 'copy',
        overwrite: params.force_download

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
        echo "[CONVERT] Streaming and converting CRAM to BAM: ${file_path}"
        samtools version
        # Stream CRAM from S3 and convert to BAM
        # samtools can read directly from S3 URIs with proper AWS credentials
        AWS_PROFILE=mdic; samtools view -b -T ${ref} -@${task.cpus} -o ${output_file_name} ${file_path}

        if [[ ! -f "${output_file_name}" ]]; then
            echo "ERROR: Failed to convert CRAM ${file_path}"
            exit 1
        fi

        echo "[CONVERT] Successfully converted to BAM"
        ls -lh "${output_file_name}"
        """
    } else if (is_s3) {
        file(file_path).copyTo(output_file_name)
    } else {
        """
        echo "[LOCAL] Using local BAM file: ${file_path}"
        ln -s "${file(file_path, checkIfExists: true)}" "${output_file_name}"
        ls -lh "${output_file_name}"
        """
    }
}

// =============================================================================
// Process: HANDLE_INDEX
// =============================================================================

process HANDLE_INDEX {
    tag "${sample_id}"

    container 'hkubal/repun:latest'

    // Publish to cache
    publishDir "${params.cache_dir}/s3_downloads",
        mode: 'copy',
        overwrite: params.force_download

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
        echo "[INDEX] Generating BAI index from converted BAM: ${bam_file}"
        samtools index -@ ${task.cpus} ${bam_file} ${output_file_name}

        if [[ ! -f "${output_file_name}" ]]; then
            echo "ERROR: Failed to generate index for ${bam_file}"
            exit 1
        fi

        ls -lh "${output_file_name}"
        """
    } else if (is_s3) {
        file(index_path).copyTo(output_file_name)
    } else {
        """
        echo "[LOCAL] Using local BAI index: ${index_path}"
        ln -s "${file(index_path, checkIfExists: true)}" "${output_file_name}"
        ls -lh "${output_file_name}"
        """
    }
}

// =============================================================================
// Workflow
// =============================================================================

workflow {

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
    DOWNLOAD_AND_CONVERT_ALIGNMENT(
        alignment_files_ch,
        ref_file,
        ref_fai_file
    )

    // Handle index files (download .bai or generate from converted BAM if source was .crai)
    // Join with alignment files so we have the BAM available for indexing
    index_files_ch
        .join(DOWNLOAD_AND_CONVERT_ALIGNMENT.out, by: 0)  // Join by sample_id
        .multiMap { sample_id, idx_path, idx_platform, idx_output, bam_file, bam_platform, bam_output ->
            index_input: tuple(sample_id, idx_path, idx_platform, idx_output)
            bam_input: tuple(sample_id, bam_file, bam_platform, bam_output)
        }
        .set { index_inputs }

    HANDLE_INDEX(
        index_inputs.index_input,
        index_inputs.bam_input
    )

    // Combine alignment and index files for RUN_REPUN
    DOWNLOAD_AND_CONVERT_ALIGNMENT.out
        .join(HANDLE_INDEX.out, by: 0)  // Join by sample_id
        .map { sample_id, bam_file, bam_platform, bam_output, bai_file, bai_platform, bai_output ->
            tuple(sample_id, bam_file, bai_file, bam_platform, bam_output)
        }
        .set { samples_prepared_ch }

    // Run repun on each sample
    RUN_REPUN(
        samples_prepared_ch,
        ref_file,
        ref_fai_file,
        truth_file
    )
}
