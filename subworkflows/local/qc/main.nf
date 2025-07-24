include { FASTP } from '../../../modules/nf-core/fastp/main'

workflow QC {
    take:
    reads // [ val(meta), path(reads) ]

    main:
    def ch_versions = Channel.empty()

    input_reads = reads.map{ meta, reads_ ->
        [meta, (reads_.size()==1) ? reads_[0] : reads_]
    }

    FASTP(
        input_reads,
        [],
        false,
        false,
        false,
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    output_reads = FASTP.out.reads.map{ meta, reads ->
        [meta, (reads instanceof Collection) ? reads : [reads]]
    }

    emit:
    fastq = output_reads
    fastp_json = FASTP.out.json
    versions = ch_versions
}
