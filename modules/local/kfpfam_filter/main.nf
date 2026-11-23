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
    val num_models

    output:
    tuple val(meta), path("${prefix}.seeds.tsv"), emit: seeds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reader = orfs.name.endsWith('.gz') ? "gunzip -c ${orfs}" : "cat ${orfs}"
    // The filter's reader takes one line per sequence and asserts it, while `seqkit translate`
    // wraps at 60 columns, so the ORFs are flattened first. HMMER_PAIRS reads the wrapped file.
    //
    // The filter writes (ORF, model index, diagonal, gapped), indices into prep's models.tsv, and
    // only pairs scoring at least FILTER_MIN_GAPPED, which is the bar HMMER then re-scores at.
    // HMMER_PAIRS wants nail's seed format, model name first, so the index becomes models.tsv's
    // NAME column: the same name pyhmmer reads from the HMM file. So prep must be built from the HMM
    // file HMMER_PAIRS reads; a different number of models means it wasn't, and is caught here,
    // before the search. HMMER_PAIRS catches a release that has as many models but other names.
    """
    models=\$((\$(wc -l < ${prep}/models.tsv) - 1))
    if [ "\$models" -ne ${num_models} ]; then
        echo "${prep} holds \$models models and the pipeline's Pfam ${num_models}: rebuild it from the same Pfam-A.hmm" >&2
        exit 1
    fi

    ${reader} \\
        | awk '/^>/ { if (NR > 1) printf "\\n"; print; next } { printf "%s", \$0 } END { printf "\\n" }' \\
        > orfs.flat.faa

    RAYON_NUM_THREADS=${task.cpus} FILTER_MIN_GAPPED=${min_gapped} ${args} \\
        filter join ${prep} ${index} ${bar} ${gate} ${block} orfs.flat.faa pairs.tsv

    awk -F'\\t' 'NR == FNR { if (FNR > 1) name[FNR - 2] = \$2; next } { print name[\$2] "\\t" \$1 }' \\
        ${prep}/models.tsv pairs.tsv > ${prefix}.seeds.tsv

    rm -f pairs.tsv orfs.flat.faa

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
