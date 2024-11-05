/*
 * Generate MultiQC report
 */
process MULTIQC {

    container 'community.wave.seqera.io/library/multiqc:1.24.1--789bc3917c8666da'
    conda "bioconda::multiqc=1.24.1"

    publishDir "$params.outdir", mode: 'copy'

    input:
        path input_files
        val cohort_name

    output:
        path "${cohort_name}_multiqc_report.html", emit: report

    script:
    """
    multiqc \\
        --force \\
        -o . \\
        -n ${cohort_name}_multiqc_report.html \\
        --clean-up \\
        ${input_files}
    """
}