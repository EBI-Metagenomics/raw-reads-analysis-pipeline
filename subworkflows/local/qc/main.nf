include { FASTP } from '../../../modules/nf-core/fastp/main'

workflow QC {
    take:
    input_reads // [ val(meta), path(reads) ]

    main:
    def ch_versions = Channel.empty()

    FASTP(
        input_reads,
        [],
        false,
        false,
        false,
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    emit:
    fastq = FASTP.out.reads
    fastp_json = FASTP.out.json
    versions = ch_versions
}
