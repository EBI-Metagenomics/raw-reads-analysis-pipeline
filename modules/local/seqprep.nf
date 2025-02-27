/*
 * SeqPrep: overlap reads
*/

process SEQPREP {

    publishDir "${params.outdir}/qc/seqprep", mode: 'copy'

    label 'seqprep'
    tag "${meta.id}"
    container 'quay.io/biocontainers/seqprep:1.3.2--hed695b0_4'

    input:
        tuple val(meta), val(decontam_reads), val(sample_d)

    output:
        tuple val(meta), path("${meta.id}_merged.fastq.gz"), path("${meta.id}_forward_unmerged.fastq.gz"), path("${meta.id}_reverse_unmerged.fastq.gz")

    script:
    def input_reads = "";
    if (decontam_reads[0].name.contains("_1")) {
        input_reads = "-f ${decontam_reads[0]} -r ${decontam_reads[1]}"
    } else {
        input_reads = "-f ${decontam_reads[1]} -r ${decontam_reads[0]}"
    }
    """
    SeqPrep \
    ${input_reads} \
    -1 ${meta.id}_forward_unmerged.fastq.gz \
    -2 ${meta.id}_reverse_unmerged.fastq.gz \
    -s ${meta.id}_merged.fastq.gz

    """
}

process SEQPREP_REPORT {

    publishDir "${params.outdir}/qc/seqprep", mode: 'copy'

    container 'quay.io/biocontainers/seqprep:1.3.2--hed695b0_4'
    tag "${meta.id}"
    label 'seqprep_report'

    input:
        tuple val(meta), val(reads)
    
    output:
        tuple val(meta), path("seqprep_output_report.txt")

    script:
    """
    zcat ${reads.forward_unmapped} | grep '@' | wc -l > seqprep_output_report.txt
    zcat ${reads.reverse_unmerged} | grep '@' | wc -l >> seqprep_output_report.txt
    zcat ${reads.merged} | grep '@' | wc -l >> seqprep_output_report.txt
    """
}
