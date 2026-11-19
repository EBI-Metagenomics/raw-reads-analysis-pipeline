process KFPFAM_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    // kfpfam's filter is not published anywhere, so the container is a local Singularity image:
    // build it with hpc/pfam-arms/kfpfam_setup.sh and pass --kfpfam_sif.
    container "${params.kfpfam_sif}"

    input:
    tuple val(meta), path(orfs)
    path prep
    path index
    val bar
    val gate
    val min_gapped
    val block

    output:
    tuple val(meta), path("${prefix}.seeds.tsv"), emit: seeds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def target = orfs.name.endsWith('.gz') ? orfs.name[0..-4] : orfs.name
    def unzip = orfs.name.endsWith('.gz') ? "gunzip -c ${orfs} > ${target}" : ''
    // The filter writes (ORF, model index, diagonal, gapped), indices into prep's models.tsv, and
    // only pairs scoring at least FILTER_MIN_GAPPED, which is the bar HMMER then re-scores at.
    // HMMER_PAIRS wants nail's seed format, model name first, so the index becomes models.tsv's
    // NAME column: the same name pyhmmer reads from the HMM file.
    """
    ${unzip}

    RAYON_NUM_THREADS=${task.cpus} FILTER_MIN_GAPPED=${min_gapped} ${args} \\
        filter join ${prep} ${index} ${bar} ${gate} ${block} ${target} pairs.tsv

    awk -F'\\t' 'NR == FNR { if (FNR > 1) name[FNR - 2] = \$2; next } { print name[\$2] "\\t" \$1 }' \\
        ${prep}/models.tsv pairs.tsv > ${prefix}.seeds.tsv

    rm -f pairs.tsv ${orfs.name.endsWith('.gz') ? target : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kfpfam: \$(cat /opt/kfpfam/VERSION 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kfpfam: stub
    END_VERSIONS
    """
}
