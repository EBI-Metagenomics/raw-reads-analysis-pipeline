process HMMER_PAIRS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyhmmer_pyarrow:2f9cd8cb7474af9e'
        : 'community.wave.seqera.io/library/python_pyhmmer_pyarrow:15ad25b8fbebe783'}"

    input:
    tuple val(meta), path(orfs), path(seeds)
    path hmm
    val num_models

    output:
    tuple val(meta), path("${prefix}.domtbl.gz"), emit: domtbl
    tuple val(meta), path("${prefix}_alignments.parquet"), emit: alignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def script = file("${moduleDir}/bin/hmmer_pairs.py")
    """
    python ${script} \\
        --orfs ${orfs} \\
        --seeds ${seeds} \\
        --hmm ${hmm} \\
        --num-models ${num_models} \\
        --threads ${task.cpus} \\
        --domtbl ${prefix}.domtbl.gz \\
        --parquet ${prefix}_alignments.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: \$(python -c 'import pyhmmer; print(pyhmmer.__version__)')
        pyarrow: \$(python -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip > ${prefix}.domtbl.gz
    touch ${prefix}_alignments.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: stub
        pyarrow: stub
    END_VERSIONS
    """
}
