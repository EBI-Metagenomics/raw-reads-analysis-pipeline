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
            .map { meta, reads_ ->
                def reads__ = reads_ instanceof Collection ? reads_[0] : reads_
                [meta] + [reads__]
            }
            .splitFastq(
                by: params.decontam_long_chunksize,
                elem: 1,
                file: true
            )
            .groupTuple()
            .mix(
                pe_se_reads.pe
                .map { meta, reads_ -> [meta] + reads_ }
                .splitFastq(
                    by: params.decontam_long_chunksize,
                    elem: [1,2],
                    file: true
                )
                .map { meta, reads1, reads2 -> tuple(meta, [reads1, reads2])}
                .groupTuple()
            )
            .flatMap {
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
                [meta, bam, !meta.single_end, ""]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

        chunked_decontaminated_reads = DECONTAMBAM.out.unmapped_reads

        concat_ch = chunked_decontaminated_reads
            .groupTuple()
            .flatMap {
                meta, reads_ ->
                def stacked_reads = reads_
                if (meta.single_end) {
                    stacked_reads = [reads_ instanceof Collection ? reads_ : [reads_]]
                } else {
                    stacked_reads = (reads_[0] instanceof Collection ? reads_ : [reads_]).transpose()
                }
                return stacked_reads.indexed().collect{
                    idx, reads_stack ->
                    tuple(groupKey(meta, stacked_reads.size()), "${meta.id}_${idx+1}.fastq.gz", reads_stack)
                }
            }
        CONCATENATE(concat_ch)

        decontaminated_reads = CONCATENATE.out.concatenated_file.groupTuple()

    }
    else {
        decontaminated_reads = reads
    }

    emit:
    decontaminated_reads = decontaminated_reads
    versions = ch_versions
}
