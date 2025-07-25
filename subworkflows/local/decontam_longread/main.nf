include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { DECONTAMBAM } from '../../../modules/local/decontambam/main'
include { CHUNKFASTX} from  '../../../modules/local/chunkfastx/main'
include { CONCATENATE} from  '../../../modules/local/concatenate/main'

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

    if (params.remove_host) {
        decontaminated_reads = input_reads.flatMap { meta, fastas ->
            def reads_n = fastas.size()
            return [fastas.indices, fastas].transpose().collect {
                idx, fasta ->
                tuple(meta + ['read_idx': idx, 'reads_n': reads_n], fasta)
            }
        }

        CHUNKFASTX(decontaminated_reads)
        chunked_reads = CHUNKFASTX.out.chunked_reads.flatMap {
            meta, chunks ->
            def chunks_ = chunks instanceof Collection ? chunks : [chunks]
            return chunks_.collect {
                chunk ->
                tuple(groupKey(meta, chunks_.size()), chunk)
            }
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

        DECONTAMBAM(
            MINIMAP2_ALIGN.out.bam.map { meta, bam ->
                [meta, bam, false, "${meta.read_idx+1}"]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

        chunked_decontaminated_reads = DECONTAMBAM.out.unmapped_reads

        CONCATENATE(
            chunked_decontaminated_reads
            .groupTuple()
            .map { meta, reads ->
                tuple(meta, "${meta.id}.${meta.read_idx+1}.fq.gz", reads)
            }
        )
        decontaminated_reads = CONCATENATE.out.concatenated_file

        decontaminated_reads = decontaminated_reads.map {
            meta, fasta ->
            def meta_ = meta - ['read_idx': meta.read_idx, 'reads_n': meta.reads_n]
            return tuple(
                groupKey(meta_, meta.reads_n),
                fasta
            )
        }
        .groupTuple()

    }
    else {
        decontaminated_reads = input_reads
    }

    emit:
    decontaminated_reads = decontaminated_reads
    versions = ch_versions
}
