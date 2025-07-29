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
            .splitFastq(
                by: params.decontam_longread_chunksize,
                pe: false,
                file: true
            )
            .mix(
                pe_se_reads.pe.splitFastq(
                    by: params.decontam_longread_chunksize,
                    pe: true,
                    file: true
                )
            )
            .flatMap {
                meta, chunks ->
                return chunks.collect {
                    chunk ->
                    tuple(groupKey(meta, chunks.size()), chunk)
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

        CONCATENATE(
            chunked_decontaminated_reads
            .groupTuple()
            .flatMap {
                meta, reads ->
                return reads.transpose()
                    .collect{ reads_ -> tuple(groupKey(meta), "${meta.id}.fastq.gz", reads_) }
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
