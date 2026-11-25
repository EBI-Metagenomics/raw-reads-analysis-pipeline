process HMMER_PAIRS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyhmmer_pyarrow:2f9cd8cb7474af9e'
        : 'community.wave.seqera.io/library/python_pyhmmer_pyarrow:15ad25b8fbebe783'}"

    input:
    tuple val(meta), path(orfs), path(hmm), path(seeds)
    val num_models
    val alignments // whether to write them
    val weak_evalue // the alignments' E-value cut for domains under the gathering thresholds too; 0 for none

    output:
    tuple val(meta), path("${prefix}.domtbl.gz"), emit: domtbl
    tuple val(meta), path("${prefix}_alignments.parquet"), emit: alignments, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def script = file("${moduleDir}/bin/hmmer_pairs.py")
    def parquet = alignments ? "--parquet ${prefix}_alignments.parquet --weak-evalue ${weak_evalue}" : ''
    // an unchunked Pfam must hold every model the seeds name
    def whole = meta.hmm_chunks == 1 ? '--whole' : ''
    """
    python ${script} \\
        --orfs ${orfs} \\
        --seeds ${seeds} \\
        --hmm ${hmm} \\
        --num-models ${num_models} ${whole} \\
        --threads ${task.cpus} \\
        --domtbl ${prefix}.domtbl.gz ${parquet}

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
    ${alignments ? "touch ${prefix}_alignments.parquet" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: stub
        pyarrow: stub
    END_VERSIONS
    """
}
