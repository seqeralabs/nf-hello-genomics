/*
 * Generate statistics with bcftools stats
 */
process BCFTOOLS_STATS {

    container 'community.wave.seqera.io/library/bcftools:1.20--a7f1d9cdda56cc93'
    conda "bioconda::bcftools=1.20"

    input:
        path vcf_file

    output:
        path "${vcf_file}.stats", emit: stats

    script:
    """
    bcftools stats ${vcf_file} > ${vcf_file}.stats
    """
}