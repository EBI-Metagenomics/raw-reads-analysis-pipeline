include { BWAMEM2_MEM as BWAMEM2_ALIGN_PHIX } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN_HOST } from '../../../modules/nf-core/bwamem2/mem/main'
include { DECONTAMBAM as DECONTAMBAM_PHIX } from '../../../modules/local/decontambam/main'
include { DECONTAMBAM as DECONTAMBAM_HOST } from '../../../modules/local/decontambam/main'

workflow DECONTAM_SHORTREAD {
    take:
    reads  // [ val(meta), path(reads) ]
    reference_genome  // [ val(meta2), path(reference_genome_index_root) ]
    phix_genome  // [ val(meta3), path(phix_index_root) ]

    main:
    def ch_versions = Channel.empty()

    phix_genome_index = phix_genome
        .map{ meta, fp ->
              [meta, files("${fp}/${meta.base_dir}/${meta.files.bwa_index_prefix}.*")] }
        .first()
    phix_genome_fasta = phix_genome
        .map{ meta, fp ->
              [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")] }
        .first()
    reference_genome_index = reference_genome
        .map{ meta, fp ->
              [meta, files("${fp}/${meta.base_dir}/${meta.files.bwa_index_prefix}.*")] }
        .first()
    reference_genome_fasta = reference_genome
        .map{ meta, fp ->
              [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")] }
        .first()

    // phix_genome_index.view{ "phix_genome_index - ${it}" }
    // phix_genome_fasta.view{ "phix_genome_fasta - ${it}" }
    // reference_genome_index.view{ "reference_genome_index - ${it}" }
    // reference_genome_fasta.view{ "reference_genome_fasta - ${it}" }

    if (params.remove_human_phix) {

        BWAMEM2_ALIGN_PHIX(
            reads,
            phix_genome_index,
            phix_genome_fasta,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_PHIX.out.versions)

        DECONTAMBAM_PHIX(
            BWAMEM2_ALIGN_PHIX.out.bam.map{ meta, bam -> [meta, bam, meta.single_end==false] }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM_PHIX.out.versions)

        decontaminated_reads = DECONTAMBAM_PHIX.out.unmapped_reads
    }
    else {
        decontaminated_reads = reads
    }

    if (reference_genome != null) {

        BWAMEM2_ALIGN_HOST(
            decontaminated_reads,
            reference_genome_index,
            reference_genome_fasta,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_HOST.out.versions)

        DECONTAMBAM_HOST(
            BWAMEM2_ALIGN_HOST.out.bam.map{ meta, bam -> [meta, bam, meta.single_end==false] }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM_HOST.out.versions)

        decontaminated_reads = DECONTAMBAM_HOST.out.unmapped_reads
    }

    emit:
    decontaminated_reads = decontaminated_reads
    versions = ch_versions
}
