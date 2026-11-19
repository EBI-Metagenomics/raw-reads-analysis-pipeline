process SEQKIT_SHUFFLE_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fasta, stageAs: 'input/*')
    val subsample_n
    val shuffle

    output:
    tuple val(meta), path("${prefix}.*"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ? "${task.ext.prefix}.${meta.id}" : "${meta.id}"
    def call_gzip = fasta.toString().endsWith('.gz') ? "| gzip -c ${args2}" : ''
    extension = fasta.toString().endsWith('.gz') ? '.gz' : ''
    // Seeded reservoir sampling: the same reads always give the same subsample, so runs and
    // arms of a comparison search the same draw. `shuf` can't do this here, since busybox's
    // has no --random-source. It also holds only the kept names, not every name.
    def seed = task.ext.seed ?: 23
    def sample_cmd = "awk -v n=${subsample_n} -v seed=${seed} '" +
        'BEGIN { srand(seed) } ' +
        '{ if (NR <= n) r[NR] = $0; else { i = int(rand() * NR) + 1; if (i <= n) r[i] = $0 } } ' +
        "END { for (j = 1; j <= n && j <= NR; j++) print r[j] }'"
    def head_cmd = shuffle ? sample_cmd : "head -n ${subsample_n}"

    """
    seqkit seq --name --only-id ${fasta} \\
    | ${head_cmd} \\
    > seq_ids

    seqkit faidx --region-file seq_ids --full-head ${fasta} \\
    ${call_gzip} \\
    > ${prefix}.fasta${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ? "${task.ext.prefix}.${meta.id}" : "${meta.id}"
    extension = fasta.toString().endsWith('.gz') ? '.gz' : ''
    """
    touch ${prefix}.fasta${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
