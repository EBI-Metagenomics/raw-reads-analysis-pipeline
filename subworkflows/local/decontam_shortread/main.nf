include { BWAMEM2_MEM_SAMTOOLS_BAM2FQ as DECONTAM_PHIX } from '../../../modules/local/bwamem2_mem_samtools_bam2fq/main'
include { BWAMEM2_MEM_SAMTOOLS_BAM2FQ as DECONTAM_HOST } from '../../../modules/local/bwamem2_mem_samtools_bam2fq/main'
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
            [meta, files("${fp}/${meta.files.bwa_index_prefix}.*")]
        }
        .first()

    host_genome_index = host_genome
        .map { meta, fp ->
            [meta, files("${fp}/${meta.files.bwa_index_prefix}.*")]
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
        
        // chunk single-end and paired-end reads seperately
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
            .groupTuple()  // group to count the number of chunks
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
            .flatMap {  // un-group with information about the number of chunks (to aid re-grouping)
                meta, chunks ->
                def chunks_ = chunks instanceof Collection ? chunks : [chunks]
                return chunks_.collect {
                    chunk ->
                    tuple(groupKey(meta, chunks_.size()), chunk)
                }
            }

        if (params.remove_phix) {
            
            decontam_phix_ch = chunked_decontaminated_reads.branch { 
                meta, _reads ->
                illumina: meta.instrument_platform=='ILLUMINA'
                not_illumina : !(meta.instrument_platform=='ILLUMINA')
            }

            DECONTAM_PHIX(
                decontam_phix_ch.illumina
                    .map{ meta, fastq -> [meta, fastq, !meta.single_end] },
                phix_genome_index,
            )
            ch_versions = ch_versions.mix(DECONTAM_PHIX.out.versions)

            chunked_decontaminated_reads = DECONTAM_PHIX.out.reads
                .mix(decontam_phix_ch.not_illumina)
        }

        if (host_genome != null) {

            DECONTAM_HOST(
                chunked_decontaminated_reads
                    .map{ meta, fastq -> [meta, fastq, !meta.single_end] },
                host_genome_index,
            )
            ch_versions = ch_versions.mix(DECONTAM_HOST.out.versions)

            chunked_decontaminated_reads = DECONTAM_HOST.out.reads
        }

        // join the chunked reads back together
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
