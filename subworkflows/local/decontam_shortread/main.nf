include { BWAMEM2_MEM as BWAMEM2_ALIGN_PHIX } from '../../../modules/nf-core/bwamem2/mem'
include { BWAMEM2_MEM as BWAMEM2_ALIGN_HOST } from '../../../modules/nf-core/bwamem2/mem'
include { SAMTOOLS_BAM2FQ as SAMTOOLS_BAM2FQ_PHIX } from '../../../modules/ebi-metagenomics/samtools/bam2fq/main'
include { SAMTOOLS_BAM2FQ as SAMTOOLS_BAM2FQ_HOST } from '../../../modules/ebi-metagenomics/samtools/bam2fq/main'

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

    phix_genome_index.view{ "phix_genome_index - ${it}" }
    phix_genome_fasta.view{ "phix_genome_fasta - ${it}" }
    reference_genome_index.view{ "reference_genome_index - ${it}" }
    reference_genome_fasta.view{ "reference_genome_fasta - ${it}" }

    if (params.remove_human_phix) {

        BWAMEM2_ALIGN_PHIX(
            reads,
            phix_genome_index,
            phix_genome_fasta,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_PHIX.out.versions)

        SAMTOOLS_BAM2FQ_PHIX(BWAMEM2_ALIGN_PHIX.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ_PHIX.out.versions)

        decontaminated_reads = SAMTOOLS_BAM2FQ_PHIX.out.reads
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

        SAMTOOLS_BAM2FQ_HOST(BWAMEM2_ALIGN_HOST.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ_HOST.out.versions)

        decontaminated_reads = SAMTOOLS_BAM2FQ_HOST.out.reads
    }

    emit:
    qc_reads = decontaminated_reads
    versions = ch_versions
}
