include { FASTP } from "${projectDir}/modules/ebi-metagenomics/fastp/main"

workflow  READSTRIM {

    take:
    ch_reads    // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    FASTP ( ch_reads, params.save_trimmed_fail, false )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    emit:
    reads               = FASTP.out.reads  // channel: [ val(meta), [ fastq ] ]
    fastp_summary_json  = FASTP.out.json   // channel: [ val(meta), [ json ] ]
    versions            = ch_versions      // channel: [ versions.yml ]
}
