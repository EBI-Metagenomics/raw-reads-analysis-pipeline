process FRAGGENESCANRS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/frag_gene_scan_rs:1.1.0--h4349ce8_0'
        : 'biocontainers/frag_gene_scan_rs:1.1.0--h4349ce8_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.out"), emit: calls
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_fgs"
    def read_cmd = reads.name.endsWith('.gz') ? "gunzip -c ${reads}" : "cat ${reads}"
    // Only coordinates are used: FragGeneScanRs 1.1.0's proteins repeat a codon at most read ends
    """
    ${read_cmd} | FragGeneScanRs \\
        ${args} \\
        -t illumina_5 \\
        -w 0 \\
        -p ${task.cpus} \\
        -m ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fraggenescanrs: \$(FragGeneScanRs --version | sed 's/FragGeneScanRs //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_fgs"
    """
    touch ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fraggenescanrs: stub
    END_VERSIONS
    """
}
