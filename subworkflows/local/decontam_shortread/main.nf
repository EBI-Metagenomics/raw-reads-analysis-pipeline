include { BWAMEM2_MEM as BWAMEM2_ALIGN_PHIX } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN_HOST } from '../../../modules/nf-core/bwamem2/mem/main'
include { DECONTAMBAM as DECONTAMBAM_PHIX } from '../../../modules/local/decontambam/main'
include { DECONTAMBAM as DECONTAMBAM_HOST } from '../../../modules/local/decontambam/main'
include { CHUNKFASTX} from  '../../../modules/local/chunkfastx/main'
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

        decontaminated_reads = reads.flatMap { meta, fastas ->
            def reads_n = fastas.size()
            return [fastas.indices, fastas].transpose().collect {
                idx, fasta ->
                tuple(meta + ['read_idx': idx, 'reads_n': reads_n], fasta)
            }
        }

        CHUNKFASTX(decontaminated_reads)
        chunked_decontaminated_reads = CHUNKFASTX.out.chunked_reads.flatMap {
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
                    [meta, bam, false, "${meta.read_idx+1}"]
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
                    [meta, bam, false, "${meta.read_idx+1}"]
                }
            )
            ch_versions = ch_versions.mix(DECONTAMBAM_HOST.out.versions)

            chunked_decontaminated_reads = DECONTAMBAM_HOST.out.unmapped_reads
        }

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
        decontaminated_reads = reads
    }


    emit:
    decontaminated_reads = decontaminated_reads
    versions = ch_versions
}
