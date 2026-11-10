process MERGEPARQUET {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyhmmer_pyarrow:2f9cd8cb7474af9e'
        : 'community.wave.seqera.io/library/python_pyhmmer_pyarrow:15ad25b8fbebe783'}"

    input:
    tuple val(meta), path(parquets, stageAs: "inputs/?/*")

    output:
    tuple val(meta), path("${prefix}_alignments.parquet"), emit: parquet
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Sorted and written as HMMER_PAIRS does, so a merged file reads like an unchunked one
    """
    python - ${parquets} <<'PY'
    import sys
    import pyarrow.parquet as pq

    files = sorted(sys.argv[1:])
    assert files, "no Parquet inputs staged"
    table = pq.read_table(files)
    table = table.sort_by([("pfam", "ascending"), ("read", "ascending")])
    pq.write_table(
        table, "${prefix}_alignments.parquet",
        compression="zstd", compression_level=9, use_dictionary=["pfam"],
    )
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyarrow: \$(python -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_alignments.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyarrow: stub
    END_VERSIONS
    """
}
