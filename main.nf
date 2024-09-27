/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'
    conda "bioconda::samtools=1.20"

    input:
        path input_bam

    output:
        tuple path(input_bam, includeInputs: true), path("${input_bam}.bai")

    """
    samtools index '$input_bam'

    """
}

/*
 * Call variants with GATK HapolotypeCaller in GVCF mode
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"
    conda "bioconda::gatk4=4.5.0.0"

    input:
        tuple path(input_bam), path(input_bam_index)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.g.vcf"
        path "${input_bam}.g.vcf.idx"

    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.g.vcf \
        -L ${interval_list} \
        -ERC GVCF
    """
}

/*
 * Consolidate GVCFs and apply joint genotyping analysis
 */
process GATK_JOINTGENOTYPING {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"
    conda "bioconda::gatk4=4.5.0.0"

    input:
        path vcfs
        path idxs
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${cohort_name}.joint.vcf"
        path "${cohort_name}.joint.vcf.idx"

    script:
    def input_vcfs = vcfs.collect { "-V ${it}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${input_vcfs} \
        --genomicsdb-workspace-path ${cohort_name}_gdb \
        -L ${interval_list}

    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${cohort_name}_gdb \
        -O ${cohort_name}.joint.vcf \
        -L ${interval_list}
    """
}

/*
 * Generate statistics with bcftools stats
 */
process BCFTOOLS_STATS {

    container 'community.wave.seqera.io/library/bcftools:1.20--a7f1d9cdda56cc93'
    conda "bioconda::bcftools=1.20"

    input:
        path vcf_file

    output:
        path "${vcf_file}.stats"

    """
    bcftools stats ${vcf_file} > ${vcf_file}.stats
    """
}

/*
 * Generate MultiQC report
 */
process MULTIQC {

    container 'community.wave.seqera.io/library/multiqc:1.24.1--789bc3917c8666da'
    conda "bioconda::multiqc=1.24.1"

    publishDir "${params.outdir}", mode: 'copy'

    input:
        path input_files
        val cohort_name

    output:
        path "${params.cohort_name}_multiqc_report.html"

    """
    multiqc \\
        --force \\
        -o . \\
        -n ${cohort_name}_multiqc_report.html \\
        --clean-up \\
        ${input_files}
    """
}

workflow {

    /*
    * Pipeline parameters
    */


    // Create input channel from samplesheet in CSV format (via CLI parameter)
    reads_ch = Channel.fromPath(params.reads_bam)
                        .splitCsv(header: true)
                        .map{ row -> [row.id, file(row.bam)] }

    // Accessory files
    params.reference = "${workflow.projectDir}/data/ref/ref.fasta"
    params.reference_index = "${workflow.projectDir}/data/ref/ref.fasta.fai"
    params.reference_dict = "${workflow.projectDir}/data/ref/ref.dict"
    params.intervals = "${workflow.projectDir}/data/ref/intervals.bed"

    // Base name for final output file
    params.cohort_name = "family_trio"

    // Output directory
    params.outdir = "results"

    // Create input channel from list of input files in plain text
    reads_ch = Channel.fromPath(params.reads_bam).splitText()
    
    // Create channels for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )

    all_vcfs = GATK_HAPLOTYPECALLER.out[0].collect()
    all_tbis = GATK_HAPLOTYPECALLER.out[1].collect()

    // Consolidate GVCFs and apply joint genotyping analysis
    GATK_JOINTGENOTYPING(
        all_vcfs,
        all_tbis,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )

    BCFTOOLS_STATS(
        GATK_JOINTGENOTYPING.out[0]
    )

    MULTIQC(
        BCFTOOLS_STATS.out.collect(),
        params.cohort_name
    )
}