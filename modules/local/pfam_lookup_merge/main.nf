process PFAM_LOOKUP_MERGE {
    tag "${meta.id}"
    label 'process_single'

    // The same image as PFAM_LOOKUP: bundle.py merge joins its parts into the profile one task would write.
    container "${params.pfam_lookup_sif}"

    input:
    tuple val(meta), path(parts, stageAs: 'parts/*')
    path bundle
    val floor

    output:
    tuple val(meta), path("${prefix}.reads.parquet"), emit: reads
    tuple val(meta), path("${prefix}.*.tsv"), emit: profiles, optional: true // one per label scheme of the bundle
    tuple val(meta), path("${prefix}.summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def script = task.ext.bundle_py ?: '/opt/lookup/bundle.py'
    prefix = task.ext.prefix ?: "${meta.id}_pfam_lookup"
    // a read's best call may be in another part than its other alignments, so it is chosen here, over all of them
    """
    python3 ${script} merge \\
        --bundle ${bundle} \\
        --parts ${parts} \\
        --out ${prefix} \\
        --floor ${floor}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polars: \$(python3 -c 'import polars; print(polars.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_pfam_lookup"
    """
    touch ${prefix}.reads.parquet ${prefix}.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polars: stub
    END_VERSIONS
    """
}
