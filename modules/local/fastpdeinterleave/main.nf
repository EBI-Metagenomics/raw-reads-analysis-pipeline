process FASTPDEINTERLEAVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/889a182b8066804f4799f3808a5813ad601381a8a0e3baa4ab8d73e739b97001/data' :
        'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def prep_in1 = "[ ! -f  ${prefix}.fastq.gz ] && ln -sf ${reads} ${prefix}.fastq.gz"
    def in1_cmd = "--in1 ${prefix}.fastq.gz"

    def out_fq1 = "--out1 ${prefix}_1.fastp.fastq.gz"
    def out_fq2 = "--out2 ${prefix}_2.fastp.fastq.gz"

    def interleaved_cmd = "--interleaved_in"

    """
    $prep_in1

    fastp \\
        $in1_cmd \\
        $out_fq1 \\
        $out_fq2 \\
        --thread $task.cpus \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        $interleaved_cmd \\
        -A -G -Q -L \\
        $args \\
        2>| >(tee ${prefix}.fastp.log >&2) \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def touch_reads         = meta.single_end ? "echo '' | gzip > ${prefix}.fastp.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastp.fastq.gz ; echo '' | gzip > ${prefix}_2.fastp.fastq.gz"
    """
    $touch_reads
    touch "${prefix}.fastp.json"
    touch "${prefix}.fastp.html"
    touch "${prefix}.fastp.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
