include { BBMAP_REFORMAT } from '../../../modules/local/bbmap/reformat/main'
include { FASTP as FASTP_SE } from '../../../modules/nf-core/fastp/main'
include { FASTP as FASTP_PE } from '../../../modules/nf-core/fastp/main'

workflow READSMERGE {
    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    // split single-end and paired-end reads
    ch_se_reads = ch_reads.filter { meta, _reads -> meta.single_end }
    ch_pe_reads = ch_reads.filter { meta, _reads -> !meta.single_end }

    FASTP_PE(ch_pe_reads, [], false, false, true)
    FASTP_SE(ch_se_reads, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP_PE.out.versions.first())

    fastp_summary_json = FASTP_SE.out.json.mix(FASTP_PE.out.json)

    // mix back with single-end reads
    ch_all_reads = ch_se_reads.mix(FASTP_PE.out.reads_merged)

    // convert to fasta
    BBMAP_REFORMAT(ch_all_reads, 'fasta')
    ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions.first())

    emit:
    reads = ch_all_reads // channel: [ val(meta), [ fastq ] ]
    fastp_summary_json = fastp_summary_json // channel: [ val(meta), [ json ] ]
    reads_fasta = BBMAP_REFORMAT.out.reformated // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
