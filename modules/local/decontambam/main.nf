process DECONTAMBAM {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biocontainers/samtools:1.21--h50ea8bc_0'
        : 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), path(inputbam), val(split), val(fn_prefix)

    output:
    tuple val(meta), path("*_mapped{_1,_2,_interleaved}.fq.gz")  , emit: mapped_reads
    tuple val(meta), path("*_unmapped{_1,_2,_interleaved}.fq.gz"), emit: unmapped_reads
    tuple val(meta), path("*_all_summary_stats.txt")             , emit: stats
    tuple val(meta), path("*_mapped_summary_stats.txt")          , emit: mapped_stats
    tuple val(meta), path("*_unmapped_summary_stats.txt")        , emit: unmapped_stats
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def all_prefix = "${prefix}_${fn_prefix}_all"
    def mapped_prefix = "${prefix}_${fn_prefix}_mapped"
    def unmapped_prefix = "${prefix}_${fn_prefix}_unmapped"
    if (split) {
        """
        #!/bin/bash
        samtools view \\
            ${args1} \\
            -f 4 \\
            -b ${inputbam} \\
        | samtools bam2fq \\
            ${args2} \\
            -@ ${task.cpus} \\
            -1 ${unmapped_prefix}_1.fq.gz \\
            -2 ${unmapped_prefix}_2.fq.gz \\
            -0 ${unmapped_prefix}_other.fq.gz \\
            -s ${unmapped_prefix}_singleton.fq.gz \\
            -

        samtools view \\
            ${args1} \\
            -F 4 \\
            -b ${inputbam} \\
        | samtools bam2fq \\
            ${args2} \\
            -@ ${task.cpus} \\
            -1 ${mapped_prefix}_1.fq.gz \\
            -2 ${mapped_prefix}_2.fq.gz \\
            -0 ${mapped_prefix}_other.fq.gz \\
            -s ${mapped_prefix}_singleton.fq.gz \\
            -

        samtools stats ${inputbam} \\
        > ${all_prefix}_summary_stats.txt

        samtools stats -F 4 ${inputbam} \\
        > ${mapped_prefix}_summary_stats.txt

        samtools stats -f 4 ${inputbam} \\
        > ${unmapped_prefix}_summary_stats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        #!/bin/bash
        samtools view \\
            ${args1} \\
            -f 4 \\
            -b ${inputbam} \\
        | samtools bam2fq \\
            ${args2} \\
            -@ ${task.cpus} \\
            - \\
        | gzip --no-name > ${unmapped_prefix}_interleaved.fq.gz

        samtools view \\
            ${args1} \\
            -F 4 \\
            -b ${inputbam} \\
        | samtools bam2fq \\
            ${args2} \\
            -@ ${task.cpus} \\
            - \\
        | gzip --no-name > ${mapped_prefix}_interleaved.fq.gz

        samtools stats ${inputbam} \\
        > ${all_prefix}_summary_stats.txt

        samtools stats -F 4 ${inputbam} \\
        > ${mapped_prefix}_summary_stats.txt

        samtools stats -f 4 ${inputbam} \\
        > ${unmapped_prefix}_summary_stats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
