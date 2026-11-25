process PFAM_LOOKUP {
    tag "${meta.id}"
    label 'process_high_memory'

    // Published nowhere, as kfpfam's image isn't: a local Singularity image holding the lookup binary and
    // bundle.py, of kmer-functional-profiler's workflows/pfam-anchor/lookup. Build it with
    // hpc/pfam-arms/lookup_setup.sh and pass --pfam_lookup_sif.
    container "${params.pfam_lookup_sif}"

    input:
    tuple val(meta), path(alignments)
    path bundle
    val floor
    val k // rows the first stage keeps
    val second // of those, rows that go to the banded score

    output:
    tuple val(meta), path("${prefix}.reads.parquet"), emit: reads
    tuple val(meta), path("${prefix}.*.tsv"), emit: profiles, optional: true // one per label scheme of the bundle
    tuple val(meta), path("${prefix}.summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def script = task.ext.bundle_py ?: '/opt/lookup/bundle.py'
    prefix = task.ext.prefix ?: "${meta.id}_pfam_lookup"
    // bundle.py stops if the alignments are to other model versions than the bundle's Pfam
    """
    python3 ${script} profile \\
        --bundle ${bundle} \\
        --parquet ${alignments} \\
        --out ${prefix} \\
        --floor ${floor} \\
        --k ${k} \\
        --second ${second} \\
        --threads ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pfam_lookup_bundle: \$(python3 -c 'import json; m = json.load(open("${bundle}/manifest.json")); print("format %s, Pfam %s, MGnify %s" % (m["format"], m["pfam"], m["mgnify"]))')
        polars: \$(python3 -c 'import polars; print(polars.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_pfam_lookup"
    """
    touch ${prefix}.reads.parquet ${prefix}.ko.tsv ${prefix}.summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pfam_lookup_bundle: stub
        polars: stub
    END_VERSIONS
    """
}
