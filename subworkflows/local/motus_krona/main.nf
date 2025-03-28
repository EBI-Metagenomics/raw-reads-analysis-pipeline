include { KRONA_KTIMPORTTEXT } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'
include { MOTUS_PROFILE } from '../../../modules/nf-core/motus/profile/main'
include { MOTUS2KRONA } from '../../../modules/local/motus2krona/main'

workflow MOTUS_KRONA {

    take:
    samples   // channel: [ val(meta), [ fastx ] ]
    motus_db  // path

    main:
    ch_versions = Channel.empty()

    MOTUS_PROFILE(
        samples,
        motus_db
    )
    ch_versions = ch_versions.mix(MOTUS_PROFILE.out.versions)
    // MOTUS_PROFILE.out.out.view { "motus.out - ${it}" }

    MOTUS2KRONA(MOTUS_PROFILE.out.out, file("${projectDir}/bin/motus2krona.py"))
    ch_versions = ch_versions.mix(MOTUS2KRONA.out.versions)

    // KRONA
    KRONA_KTIMPORTTEXT(MOTUS2KRONA.out.krona)
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}

