include { BWAMEM2DECONTNOBAMS as HUMAN_PHIX_DECONTAMINATION } from '../../../modules/ebi-metagenomics/bwamem2decontnobams/main'
include { BWAMEM2DECONTNOBAMS as HOST_DECONTAMINATION } from '../../../modules/ebi-metagenomics/bwamem2decontnobams/main'

workflow SHORT_READS_QC {
    take:
    reads // [ val(meta), path(reads) ]
    reference_genome // [ val(meta2), path(reference_genome_index_root) ] | meta2 contains the name of the reference genome
    phix_genome // [ val(meta3), path(phix_index_root) ]

    main:
    def ch_versions = Channel.empty()

    if (params.remove_human_phix) {

        HUMAN_PHIX_DECONTAMINATION(
            reads,
            phix_genome,
        )

        ch_versions = ch_versions.mix(HUMAN_PHIX_DECONTAMINATION.out.versions)

        decontaminated_reads = HUMAN_PHIX_DECONTAMINATION.out.decont_reads
    }
    else {
        decontaminated_reads = reads
    }

    if (reference_genome != null) {

        HOST_DECONTAMINATION(
            HUMAN_PHIX_DECONTAMINATION.out.decont_reads,
            reference_genome,
        )

        ch_versions = ch_versions.mix(HOST_DECONTAMINATION.out.versions)

        decontaminated_reads = HOST_DECONTAMINATION.out.decont_reads
    }

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
