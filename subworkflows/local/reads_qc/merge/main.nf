include { FASTP } from '../../../../modules/nf-core/fastp/main'
include { SEQTK_SEQ } from '../../../../modules/ebi-metagenomics/seqtk/seq/main'

workflow READSMERGE {
    take:
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:
    ch_versions = Channel.empty()

    ch_se_reads = ch_reads.filter { meta, _reads -> meta.single_end }
    ch_pe_reads = ch_reads.filter { meta, _reads -> !meta.single_end }

    FASTP(ch_pe_reads, [], false, params.save_trimmed_fail, true)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_all_reads = ch_se_reads.mix(FASTP.out.reads_merged)

    SEQTK_SEQ(ch_all_reads)
    ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions.first())

    emit:
    reads = ch_all_reads // channel: [ val(meta), [ fastq ] ]
    fastp_summary_json = FASTP.out.json // channel: [ val(meta), [ json ] ]
    reads_fasta = SEQTK_SEQ.out.fastx // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
