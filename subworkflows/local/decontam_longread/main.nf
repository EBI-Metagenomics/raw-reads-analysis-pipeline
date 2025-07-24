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

    chunked_reads = input_reads.flatMap { meta, fasta ->
        def reads_n = fasta.size()
        [fasta.indices, fasta].transpose().collect { idx, fasta_ ->
            def chunks = fasta_.splitFasta(
                file: true,
                size: params.decontam_long_host_chunksize
            )
            return chunks.collect{chunk -> tuple(idx, chunk, chunks.size())}
        }
        .collectMany { it }
        .collect{ idx, chunk, chunksize -> tuple(groupKey(meta + ['read_idx': idx, 'reads_n': reads_n], chunksize), chunk) }
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
            [meta, bam, false, "long_read_host"]
        }
    )
    ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

    decontaminated_reads = DECONTAMBAM.out.unmapped_reads.map { meta, reads ->
        [meta + ['decontam_host_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    decontaminated_reads = decontaminated_reads.filter { meta, _reads ->
        meta.decontam_host_read_count > 0
    }

    emit:
    decontaminated_reads = decontaminated_reads
    stats = DECONTAMBAM.out.unmapped_stats
    versions = ch_versions
}
