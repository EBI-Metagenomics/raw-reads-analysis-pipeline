include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'

workflow DECONTAM_LONGREAD {
    take:
    input_reads // [ val(meta), path(reads) ]
    reference_genome // [ val(meta2), path(reference_genome) ]

    main:
    def ch_versions = Channel.empty()

    MINIMAP2_ALIGN(
        input_reads,
        reference_genome,
        "human",
        "fastq",
        true,
        "bai",
        false,
        true,
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_HUMAN.out.versions)

    decontaminated_reads = MINIMAP2_ALIGN_HUMAN.out.filtered_output

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
