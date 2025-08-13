process BWAMEM2_MEM_SAMTOOLS_BAM2FQ {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*'), val(split)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*{_1,_2,_interleaved}.fq.gz"),         emit: reads
    tuple val(meta), path("*_singleton.fq.gz"),   optional: true, emit: singleton_reads
    tuple val(meta), path("*_other.fq.gz"),       optional: true, emit: other_reads
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def database = task.ext.database ?: meta2.id
    def sam_flag = meta.single_end ? 4 : 12

    def bam2fq_cmd = """ \\
    samtools bam2fq \\
        $args3 \\
        -@ $task.cpus \\
        | gzip --no-name > ${prefix}_interleaved.fq.gz \\
    """
    if (split) {
        bam2fq_cmd = """ \\
        samtools bam2fq \\
            $args3 \\
            -@ $task.cpus \\
            -1 ${prefix}_1.fq.gz \\
            -2 ${prefix}_2.fq.gz \\
            -0 ${prefix}_other.fq.gz \\
            -s ${prefix}_singleton.fq.gz \\
        """
    }
    
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem2 mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
    | samtools view -@ ${task.cpus} -f ${sam_flag} $args2 - \\
    | ${bam2fq_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
