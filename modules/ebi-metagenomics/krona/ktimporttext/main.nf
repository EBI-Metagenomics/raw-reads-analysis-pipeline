// This krona/ktimporttext module is copied from the already-existing nf-core module (https://nf-co.re/modules/krona_ktimporttext, https://github.com/nf-core/modules/commit/8fc1d24c710ebe1d5de0f2447ec9439fd3d9d66a)
// This is because there are not currently any nf-core ways of adding modules from more than one nf-core repo
// Two slight changes compared to the original is I've added a stub to this module, and it outputs a meaningful message to HTML if the data are empty

process KRONA_KTIMPORTTEXT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1':
        'biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ('*.html'), emit: html
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if grep -qvE '^[[:space:]]*(#|$)' "$report"; then
        ktImportText \\
            $args \\
            -o "${prefix}.html" \\
            "$report"
    else
        cat > "${prefix}.html" <<EOF
<html>
<body>
<h2>No taxonomic classifications were produced for this sample.</h2>
</body>
</html>
EOF
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
    END_VERSIONS
    """
}
