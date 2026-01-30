
bcftools view -s CCDC6-RET--dragenv4.3.6 final_consensus_raw.vcf.gz | 
    bcftools annotate --set-id '%VKX' |
    bcftools annotate -x "FORMAT/AD,FORMAT/AF,FORMAT/AFDP,FORMAT/ALTHC,FORMAT/ALT_F1R2,FORMAT/ALT_F2R1,FORMAT/BaseQRankSumPS,FORMAT/ClippingRankSumPS,FORMAT/DPHC,FORMAT/FOXOG,FORMAT/MQRankSumPS,FORMAT/NBQPS,FORMAT/QSS,FORMAT/REF_F1R2,FORMAT/REF_F2R1,FORMAT/ReadPosEndDistPS,FORMAT/ReadPosRankSumPS,FORMAT/GQ,FORMAT/DP,FORMAT/PL,FORMAT/PS,FORMAT/SQ,FORMAT/F1R2,FORMAT/F2R1,FORMAT/SB,FORMAT/MB" |
    bcftools reheader --samples-list final_consensus_raw |
    bcftools view --write-index -Oz -o single_trim_format.vcf.gz
