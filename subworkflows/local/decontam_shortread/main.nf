include { BWAMEM2_MEM as BWAMEM2_ALIGN_PHIX } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN_HOST } from '../../../modules/nf-core/bwamem2/mem/main'
include { DECONTAMBAM as DECONTAMBAM_PHIX } from '../../../modules/local/decontambam/main'
include { DECONTAMBAM as DECONTAMBAM_HOST } from '../../../modules/local/decontambam/main'
include { COMBINEBAM as COMBINEBAM_PHIX } from '../../../modules/local/combinebam/main'
include { COMBINEBAM as COMBINEBAM_HOST } from '../../../modules/local/combinebam/main'

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

    decontaminated_reads = reads.flatMap { meta, fastas ->
        def reads_n = fastas.size()
        return [fastas.indices, fastas].transpose().collect {
            idx, fasta ->
            tuple(meta + ['read_idx': idx, 'reads_n': reads_n], fasta)
        }
    }

    if (params.remove_phix) {
        chunked_reads = decontaminated_reads.flatMap { meta, fasta ->
            def chunks = fasta.splitFasta(
                file: true,
                size: params.decontam_short_phix_chunksize
            )
            return chunks.collect { chunk -> tuple(groupKey(meta, chunks.size), chunk) }
        }

        chunked_reads.view{ "chunked_reads - ${it}" }

        BWAMEM2_ALIGN_PHIX(
            chunked_reads,
            phix_genome_index,
            phix_genome_fasta,
            false,
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_PHIX.out.versions)

        COMBINEBAM_PHIX(
            BWAMEM2_ALIGN_PHIX.out.bam.groupTuple()
        )

        DECONTAMBAM_PHIX(
            COMBINEBAM_PHIX.out.concatenated_result.map { meta, bam ->
                [meta, bam, false, "short_read_phix_${meta.read_idx+1}"]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM_PHIX.out.versions)


        decontaminated_reads = DECONTAMBAM_PHIX.out.unmapped_reads
        phix_stats = DECONTAMBAM_PHIX.out.unmapped_stats
    }
    else {
        phix_stats = Channel.empty()
    }

    decontaminated_reads = decontaminated_reads.map { meta, reads_ ->
        [meta + ['decontam_phix_read_count': reads_.countFastq()], reads_]
    }
    decontaminated_reads = decontaminated_reads.filter { meta, _reads ->
        meta.decontam_phix_read_count > 0
    }

    if (host_genome != null) {
        chunked_decontaminated_reads = decontaminated_reads.flatMap { meta, fasta ->
            def chunks = fasta.splitFasta(
                file: true,
                size: params.decontam_short_host_chunksize
            )
            return chunks.collect { chunk -> tuple(groupKey(meta, chunks.size), chunk) }
        }

        BWAMEM2_ALIGN_HOST(
            chunked_decontaminated_reads,
            host_genome_index,
            host_genome_fasta,
            false,
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_HOST.out.versions)

        COMBINEBAM_HOST(
            BWAMEM2_ALIGN_HOST.out.bam.groupTuple()
        )

        DECONTAMBAM_HOST(
            COMBINEBAM_HOST.out.concatenated_result.map { meta, bam ->
                [meta, bam, false, "short_read_host_${meta.read_idx+1}"]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM_HOST.out.versions)

        decontaminated_reads = DECONTAMBAM_HOST.out.unmapped_reads
        host_stats = DECONTAMBAM_HOST.out.unmapped_stats
    }
    else {
        host_stats = Channel.empty()
    }

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
    host_stats = host_stats
    phix_stats = phix_stats
    versions = ch_versions
}
