#!/usr/bin/sh
WRKDIR=/wrk/mdic_repun
REFDIR=/resources
FA=${REFDIR}/dragen_hg38_multigenome_v4.fasta \
FAIDX={REFDIR}/dragen_hg38_multigenome_v4.fasta.fai
CRAM=${WRKDIR}/data/CCDC6-RET_WGS_120x_B1_dragen_hg38_multigenome_v4_tumor.cram

vmtouch -ld ${CRAM}.crai ${FA} ${FAIDX}

export REF_PATH=/resources/ref_cache/%2s/%2s/%s
export REF_CACHE=/resources/ref_cache/%2s/%2s/%s

time docker run -it \
  -v ${WRKDIR}:${WRKDIR} \
  -v ${REFDIR}:${REFDIR} \
  repun:v0.1.4 \
  /opt/bin/repun --bam_fn ${CRAM} --ref_fn ${FA} \
	--truth_vcf_fn ${WRKDIR}/data/CCDC6-RET--dragenv4.3.6.vcf.gz \
	-t 24 --platform ilmn \
	--somatic_mode --max_af_for_somatic_unification 0.01 --vaf_threshold_for_pass 0.01 \
	--output_dir ${WRKDIR}/CCDC6-RET_WGS_120x-cram-wtouch-dragref

