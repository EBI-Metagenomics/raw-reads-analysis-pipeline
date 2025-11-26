process MINIMAP2_ALIGN_SAMTOOLS_BAM2FQ {
    tag "$meta.id"
    label 'process_low'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66dc96eff11ab80dfd5c044e9b3425f52d818847b9c074794cf0c02bfa781661/data' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*'), val(split)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*{_1,_2,_interleaved}.fq.gz"),         emit: reads
    tuple val(meta), path("*_singleton.fq.gz"),   optional: true, emit: singleton_reads
    tuple val(meta), path("*_other.fq.gz"),       optional: true, emit: other_reads
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools reset --threads ${task.cpus-1} $args3 $reads | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)
    def sam_flag = meta.single_end ? 4 : 12

    def bam2fq_cmd = """ \\
    samtools bam2fq \\
        $args3 \\
        -@ $task.cpus \\
        | gzip --no-name > ${prefix}_interleaved.fq.gz
    """
    if (split) {
        bam2fq_cmd = """ \\
        samtools bam2fq \\
            $args3 \\
            -@ $task.cpus \\
            -1 ${prefix}_1.fq.gz \\
            -2 ${prefix}_2.fq.gz \\
            -0 ${prefix}_other.fq.gz \\
            -s ${prefix}_singleton.fq.gz\\
        """
    }

    """
    $samtools_reset_fastq \\
    minimap2 \\
        $args \\
        -t $task.cpus \\
        -a \\
        $target \\
        $query \\
    | samtools view -@ ${task.cpus} -f ${sam_flag} $args2 - \\
    | ${bam2fq_cmd}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
