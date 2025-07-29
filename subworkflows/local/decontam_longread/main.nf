include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { DECONTAMBAM } from '../../../modules/local/decontambam/main'
include { CONCATENATE} from  '../../../modules/local/concatenate/main'

workflow DECONTAM_LONGREAD {
    take:
    reads // [ val(meta), path(reads) ]
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
        reads = reads.map { meta, fastqs ->
            return tuple(meta + ['reads_n': fastqs.size()], fastqs)
        }

        pe_se_reads = reads.branch{
            meta, _reads ->
            se: meta.single_end
            pe: !meta.single_end
        }
        chunked_reads = pe_se_reads.se
            .map { meta, reads_ -> [meta] + [reads_[0]] }
            .splitFastq(
                by: params.decontam_long_chunksize,
                elem: 1,
                file: true
            )
            .flatMap {
                meta, chunks ->
                def chunks_ = chunks instanceof Collection ? chunks : [chunks]
                return chunks_.collect {
                    chunk ->
                    tuple(groupKey(meta, chunks_.size()), chunk)
                }
            }
            .mix(
                pe_se_reads.pe
                .map { meta, reads_ -> [meta] + reads_ }
                .splitFastq(
                    by: params.decontam_long_chunksize,
                    elem: [1,2],
                    file: true
                )
                .flatMap {
                    meta, chunks1, chunks2 ->
                    def chunks1_ = chunks1 instanceof Collection ? chunks1 : [chunks1]
                    def chunks2_ = chunks2 instanceof Collection ? chunks2 : [chunks2]
                    return [chunks1_, chunks2_].transpose().collect {
                        chunk1, chunk2 ->
                        tuple(groupKey(meta, chunks1_.size()), [chunk1, chunk2])
                    }
                }
            )

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
                [meta, bam, !meta.single_end, ""]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

        chunked_decontaminated_reads = DECONTAMBAM.out.unmapped_reads

        CONCATENATE(
            chunked_decontaminated_reads
            .groupTuple()
            .flatMap {
                meta, reads_ ->
                def reads__ = reads_
                if (meta.single_end) {
                    reads__ = [reads_ instanceof Collection ? reads_ : [reads_]]
                } else {
                    reads__ = reads_[0] instanceof Collection ? reads_ : [reads_]
                }
                def reads_t = reads__.transpose()
                return reads_t.indexed().collect{
                    idx, reads_t_ ->
                    tuple(groupKey(meta, reads_t.size()), "${meta.id}_${idx}.fastq.gz", reads_t_)
                }
            }
        )
        decontaminated_reads = CONCATENATE.out.concatenated_file.groupTuple()

    }
    else {
        decontaminated_reads = reads
    }

    emit:
    decontaminated_reads = decontaminated_reads
    versions = ch_versions
}
