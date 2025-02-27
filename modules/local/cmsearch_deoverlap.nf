/*
 * cmsearch deoverlap Perl script
*/

process CMSEARCH_DEOVERLAP {

    label 'cmsearch_deoverlap'
    tag "${meta.id}"
    container 'quay.io/biocontainers/perl:5.22.2.1'

    input:
    tuple val(meta), path(cmsearch_matches)
    tuple val(meta_db), path(clan_information)

    output:
    tuple val(meta), path("${cmsearch_matches}.deoverlapped") 

    script:
    """
    cmsearch-deoverlap.pl --clanin ${clan_information} ${cmsearch_matches}
    """
}
