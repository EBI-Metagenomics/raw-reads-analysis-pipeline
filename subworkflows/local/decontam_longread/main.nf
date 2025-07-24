include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { DECONTAMBAM } from '../../../modules/local/decontambam/main'
include { COMBINEBAM } from '../../../modules/local/combinebam/main'
include { CHUNKFASTX } from '../../../modules/local/chunkfastx/main'
include { GZIPALL } from '../../../modules/local/gzipall/main'

workflow DECONTAM_LONGREAD {
    take:
    input_reads // [ val(meta), path(reads) ]
    reference_genome // [ val(meta2), path(reference_genome) ]

    main:
    def ch_versions = Channel.empty()

    reference_genome_index = reference_genome
        .map { meta, fp ->
            [
                meta,
                file("${fp}/${meta.base_dir}/${meta.files.index}"),
            ]
        }
        .first()

    decontaminated_reads = input_reads.flatMap { meta, fastas ->
        def reads_n = fastas.size()
        return [fastas.indices, fastas].transpose().collect {
            idx, fasta ->
            tuple(meta + ['read_idx': idx, 'reads_n': reads_n], fasta)
        }
    }

    chunked_reads = decontaminated_reads.flatMap { meta, fasta ->
        def chunks = fasta.splitFasta(
            file: true,
            size: params.decontam_long_host_chunksize
        )
        return chunks.collect { chunk -> tuple(groupKey(meta, chunks.size), chunk) }
    }

    MINIMAP2_ALIGN(
        chunked_reads,
        reference_genome_index,
        true,
        "bai",
        false,
        false,
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    COMBINEBAM(
        MINIMAP2_ALIGN.out.bam.groupTuple()
    )

    DECONTAMBAM(
        COMBINEBAM.out.concatenated_result.map { meta, bam ->
            [meta, bam, false, "long_read_host_${meta.read_idx+1}"]
        }
    )
    ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

    decontaminated_reads = decontaminated_reads.map {
        meta, fasta ->
        def meta_ = meta - ['read_idx': meta.read_idx, 'reads_n': meta.reads_n]
        return tuple(
            groupKey(meta_, meta.reads_n),
            fasta
        )
    }
    .groupTuple()


    emit:
    decontaminated_reads = decontaminated_reads
    stats = DECONTAMBAM.out.unmapped_stats
    versions = ch_versions
}
