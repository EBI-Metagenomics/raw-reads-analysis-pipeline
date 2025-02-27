
process FASTP {

    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: "*html"
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: "*json"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*fastq*"

    container 'quay.io/biocontainers/fastp:0.23.1--h79da9fb_0'
    label 'fastp'
    tag "$meta.id"

    input:
    tuple val(meta), val(fps)
    tuple val(merged_meta), val(merged_reads)

    output:
    tuple val(meta), path("${meta.id}_fastp*.fastq.gz"), path("*_fastp.json"), path("*_fastp.html"), emit: reads
    tuple val(merged_meta), path("*_merged*"), optional: true, emit: merged

    script:
    /* Handle the input reads */
    def input_reads = "";
    def output_reads = "";
    def report_name = "qc";

    if ( fps.mode == "single" ) {
        input_reads = "--in1 ${fps.reads}";
        output_reads = "--out1 ${meta.id}_fastp.fastq.gz";
    }

    if ( fps.mode == "paired" ) {
        input_reads = "--in1 ${fps.reads[0]} --in2 ${fps.reads[1]} --detect_adapter_for_pe";
        output_reads = "--out1 ${meta.id}_fastp_1.fastq.gz --out2 ${meta.id}_fastp_2.fastq.gz";
    }

    /* Optional parameters */
    def args = ""
    if ( merged_reads ) {
        args += " -m --merged_out ${meta.id}_${merged_reads}" +
        " --unpaired1 ${meta.id}.unpaired_1.fastq.gz " +
        " --unpaired2 ${meta.id}.unpaired_2.fastq.gz"
        report_name = "overlap"
    }
    args += params.fastp_params.length_filter ? " -l ${params.fastp_params.length_filter}" : "";
    args += params.fastp_params.polya_trim ? " -x ${params.fastp_params.polya_trim}" : "";
    args += params.fastp_params.qualified_quality_phred ? " -q ${params.fastp_params.qualified_quality_phred}" : "";
    args += params.fastp_params.unqualified_percent_limit ? " -u ${params.fastp_params.unqualified_percent_limit}" : "";

    """
    fastp -w ${task.cpus} \
    ${input_reads} \
    ${output_reads} \
    --json ${meta.id}_${report_name}_fastp.json \
    --html ${meta.id}_${report_name}_fastp.html \
    ${args}
    """
}
