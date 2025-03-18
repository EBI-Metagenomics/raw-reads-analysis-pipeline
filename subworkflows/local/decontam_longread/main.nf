include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { DECONTAMBAM } from '../../../modules/local/decontambam/main'

workflow DECONTAM_LONGREAD {
    take:
    input_reads // [ val(meta), path(reads) ]
    reference_genome // [ val(meta2), path(reference_genome) ]

    main:
    def ch_versions = Channel.empty()

    reference_genome_index = reference_genome
        .map{ meta, fp -> [
            meta,
            file("${fp}/${meta.base_dir}/${meta.files.index}")
        ] }
        .first()
    // reference_genome_index.view{ "minimap2_reference - ${it}" }
    // reference_genome_fasta = reference_genome
    //     .map{ meta, fp ->
    //           [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")] }
    //     .first()

    MINIMAP2_ALIGN(
        input_reads,
        reference_genome_index,
        true,
        true,
        "bai",
        false,
        false,
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // DECONTAMBAM(
    //     MINIMAP2_ALIGN.out.bam.map{ meta, bam ->
    //                                 [meta, bam, meta.single_end==false] }
    // )
    // ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

    emit:
    decontaminated_reads = MINIMAP2_ALIGN.out.unmapped_reads
    versions = ch_versions
}
