process KFPFAM_HMMER {
    tag "${meta.id}"
    label 'process_medium'

    // kfpfam's filter is not published anywhere, so the container is a local Singularity image, with
    // pyhmmer and pyarrow beside the filter: build it with hpc/pfam-arms/kfpfam_setup.sh and pass --kfpfam_sif.
    container "${params.kfpfam_sif}"

    input:
    tuple val(meta), path(orfs)
    path hmms // every chunk of a chunked Pfam: the index covers the whole database
    path prep
    path index
    val bar
    val gate
    val min_gapped
    val block
    val num_models
    val alignments // whether to write them
    val weak_evalue // the alignments' E-value cut for domains under the gathering thresholds too; 0 for none
    val forward_slack // bits under a model's threshold by Forward past which an ORF gets no domain definition; '' for none

    output:
    tuple val(meta), path("${prefix}.domtbl.gz"), emit: domtbl
    tuple val(meta), path("${prefix}_alignments.parquet"), emit: alignments, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def script = file("${moduleDir}/../hmmer_pairs/bin/hmmer_pairs.py")
    def reader = orfs.name.endsWith('.gz') ? "gunzip -c ${orfs}" : "cat ${orfs}"
    def parquet = alignments ? "--parquet ${prefix}_alignments.parquet --weak-evalue ${weak_evalue}" : ''
    def slack = forward_slack != null && forward_slack.toString() != '' ? "--forward-slack ${forward_slack}" : ''
    // One job, because the candidate pairs are far larger than anything kept: 11 GB of seeds per 4 M
    // reads at kfpfam's most sensitive setting. The filter writes each block's pairs to the pipe as
    // numbers, 8 bytes each, and HMMER searches them while the next block is filtered, so nothing is
    // written to disk and neither side holds the sample. Both use all the task's cores; the filter
    // waits when HMMER is behind.
    //
    // The filter's reader takes one line per sequence and asserts it, while `seqkit translate`
    // wraps at 60 columns, so the ORFs are flattened on the way in. hmmer_pairs.py reads the file itself.
    //
    // The filter's model numbers index prep's models.tsv, so prep must be built from the HMMs searched here:
    // hmmer_pairs.py stops if their names differ, and a different number of models is caught before the search.
    """
    set -o pipefail

    # top's prep holds a gapped bar per model (kfpfam_setup.sh); a prep without one uses min_gapped alone
    if [ -s ${prep}/gapped.txt ]; then
        export FILTER_GAPPED_FILE=${prep}/gapped.txt
    fi

    models=\$((\$(wc -l < ${prep}/models.tsv) - 1))
    if [ "\$models" -ne ${num_models} ]; then
        echo "${prep} holds \$models models and the pipeline's Pfam ${num_models}: rebuild it from the same Pfam-A.hmm" >&2
        exit 1
    fi

    ${reader} \\
        | awk '/^>/ { if (NR > 1) printf "\\n"; print; next } { printf "%s", \$0 } END { printf "\\n" }' \\
        | RAYON_NUM_THREADS=${task.cpus} FILTER_MIN_GAPPED=${min_gapped} ${args} FILTER_BINARY=1 \\
            filter join ${prep} ${index} ${bar} ${gate} ${block} /dev/stdin /dev/stdout \\
        | python3 ${script} \\
            --orfs ${orfs} \\
            --stream ${prep}/models.tsv \\
            --hmm ${hmms} \\
            --num-models ${num_models} \\
            --threads ${task.cpus} \\
            --domtbl ${prefix}.domtbl.gz ${parquet} ${slack}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kfpfam: \$(cat /opt/kfpfam/VERSION 2>/dev/null || echo unknown)
        pyhmmer: \$(python3 -c 'import pyhmmer; print(pyhmmer.__version__)')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip > ${prefix}.domtbl.gz
    ${alignments ? "touch ${prefix}_alignments.parquet" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kfpfam: stub
        pyhmmer: stub
        pyarrow: stub
    END_VERSIONS
    """
}
