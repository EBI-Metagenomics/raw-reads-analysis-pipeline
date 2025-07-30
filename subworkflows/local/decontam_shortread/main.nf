include { BWAMEM2_MEM as BWAMEM2_ALIGN_PHIX } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN_HOST } from '../../../modules/nf-core/bwamem2/mem/main'
include { DECONTAMBAM as DECONTAMBAM_PHIX } from '../../../modules/local/decontambam/main'
include { DECONTAMBAM as DECONTAMBAM_HOST } from '../../../modules/local/decontambam/main'
include { CONCATENATE} from  '../../../modules/local/concatenate/main'

workflow DECONTAM_SHORTREAD {
    take:
    reads // [ val(meta), path(reads) ]
    host_genome // [ val(meta2), path(reference_genome_index_root) ]
    phix_genome // [ val(meta3), path(phix_index_root) ]

    main:
    def ch_versions = Channel.empty()

    phix_genome_index = phix_genome
        .map { meta, fp ->
            [meta, files("${fp}/${meta.base_dir}/${meta.files.bwa_index_prefix}.*")]
        }
        .first()
    phix_genome_fasta = phix_genome
        .map { meta, fp ->
            [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")]
        }
        .first()
    host_genome_index = host_genome
        .map { meta, fp ->
            [meta, files("${fp}/${meta.base_dir}/${meta.files.bwa_index_prefix}.*")]
        }
        .first()
    host_genome_fasta = host_genome
        .map { meta, fp ->
            [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")]
       }
        .first()


    if (params.remove_phix || params.remove_host) {

        reads = reads.map { meta, fastqs ->
            return tuple(meta + ['reads_n': fastqs.size()], fastqs)
        }

        pe_se_reads = reads.branch{
            meta, _reads ->
            se: meta.single_end
            pe: !meta.single_end
        }
        chunked_decontaminated_reads = pe_se_reads.se
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

        if (params.remove_phix) {

            BWAMEM2_ALIGN_PHIX(
                chunked_decontaminated_reads,
                phix_genome_index,
                phix_genome_fasta,
                false,
            )
            ch_versions = ch_versions.mix(BWAMEM2_ALIGN_PHIX.out.versions)

            DECONTAMBAM_PHIX(
                BWAMEM2_ALIGN_PHIX.out.bam.map { meta, bam ->
                    [meta, bam, !meta.single_end, ""]
                }
            )
            ch_versions = ch_versions.mix(DECONTAMBAM_PHIX.out.versions)

            chunked_decontaminated_reads = DECONTAMBAM_PHIX.out.unmapped_reads
        }

        if (host_genome != null) {

            BWAMEM2_ALIGN_HOST(
                chunked_decontaminated_reads,
                host_genome_index,
                host_genome_fasta,
                false,
            )
            ch_versions = ch_versions.mix(BWAMEM2_ALIGN_HOST.out.versions)

            DECONTAMBAM_HOST(
                BWAMEM2_ALIGN_HOST.out.bam.map { meta, bam ->
                    [meta, bam, !meta.single_end, ""]
                }
            )
            ch_versions = ch_versions.mix(DECONTAMBAM_HOST.out.versions)

            chunked_decontaminated_reads = DECONTAMBAM_HOST.out.unmapped_reads
        }

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
