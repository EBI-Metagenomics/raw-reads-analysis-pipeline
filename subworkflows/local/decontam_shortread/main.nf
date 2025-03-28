include { BWAMEM2_MEM as BWAMEM2_ALIGN_PHIX } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN_HOST } from '../../../modules/nf-core/bwamem2/mem/main'
include { DECONTAMBAM as DECONTAMBAM_PHIX } from '../../../modules/local/decontambam/main'
include { DECONTAMBAM as DECONTAMBAM_HOST } from '../../../modules/local/decontambam/main'
include { COMBINEBAM as COMBINEBAM_PHIX } from '../../../modules/local/combinebam/main'
include { COMBINEBAM as COMBINEBAM_HOST } from '../../../modules/local/combinebam/main'
// include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_PHIX } from '../../../modules/nf-core/seqkit/split2/main'
// include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_HOST } from '../../../modules/nf-core/seqkit/split2/main'
include { CHUNKFASTX as CHUNKFASTX_PHIX } from '../../../modules/local/chunkfastx/main'
include { CHUNKFASTX as CHUNKFASTX_HOST } from '../../../modules/local/chunkfastx/main'

workflow DECONTAM_SHORTREAD {
    take:
    reads  // [ val(meta), path(reads) ]
    host_genome  // [ val(meta2), path(reference_genome_index_root) ]
    phix_genome  // [ val(meta3), path(phix_index_root) ]

    main:
    def ch_versions = Channel.empty()

    phix_genome_index = phix_genome
        .map{ meta, fp ->
              [meta, files("${fp}/${meta.base_dir}/${meta.files.bwa_index_prefix}.*")] }
        .first()
    phix_genome_fasta = phix_genome
        .map{ meta, fp ->
              [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")] }
        .first()
    host_genome_index = host_genome
        .map{ meta, fp ->
              [meta, files("${fp}/${meta.base_dir}/${meta.files.bwa_index_prefix}.*")] }
        .first()
    host_genome_fasta = host_genome
        .map{ meta, fp ->
              [meta, file("${fp}/${meta.base_dir}/${meta.files.genome}")] }
        .first()

    // phix_genome_index.view{ "phix_genome_index - ${it}" }
    // phix_genome_fasta.view{ "phix_genome_fasta - ${it}" }
    // reference_genome_index.view{ "reference_genome_index - ${it}" }
    // reference_genome_fasta.view{ "reference_genome_fasta - ${it}" }

    if (params.remove_phix) {

        // chunked_reads = reads
        //     .flatMap{ meta, fastas ->
        //         def chunks = []
        //         if(meta.single_end) {
        //             chunks = fastas.splitFasta(file: true, compress: true, size: 500.KB)
        //         } else {
        //             chunks = [
        //                 fastas[0].splitFasta(file: true, compress: true, size: 500.KB),
        //                 fastas[1].splitFasta(file: true, compress: true, size: 500.KB)
        //             ].transpose()
        //         }
        //         chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
        //     }
        // chunked_reads.view{ "chunked_reads - ${it}" }

        // SEQKIT_SPLIT2_PHIX(reads)
        CHUNKFASTX_PHIX(reads, file("${projectDir}/bin/chunk_fastx.py"))
        chunked_reads = CHUNKFASTX_PHIX.out.reads
            .flatMap{ meta, chunks ->
                if (chunks instanceof List) {
                    def grouped_chunks = [:]
                    chunks.collect{
                        chunk ->
                        def part = (chunk.name =~ /\.(chunk-\d+)\./)[0][1]
                        if (grouped_chunks.containsKey(part)) {
                            grouped_chunks[part].add(chunk)
                        } else {
                            grouped_chunks[part] = [chunk]
                        }
                    }
                    return grouped_chunks.collect{ _part, chunk -> tuple(groupKey(meta, grouped_chunks.size()), tuple(*chunk)) }
                } else {
                    return [tuple(groupKey(meta, 1), chunks)]
                }
            }
        // chunked_reads.view{ "chunked_reads - ${it}" }

        BWAMEM2_ALIGN_PHIX(
            chunked_reads,
            phix_genome_index,
            phix_genome_fasta,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_PHIX.out.versions)

        COMBINEBAM_PHIX(
            BWAMEM2_ALIGN_PHIX.out.bam.groupTuple()
        )

        DECONTAMBAM_PHIX(
            COMBINEBAM_PHIX.out.concatenated_result.map{ meta, bam ->
                [meta, bam, meta.single_end==false, "short_read_phix"]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM_PHIX.out.versions)


        decontaminated_reads = DECONTAMBAM_PHIX.out.unmapped_reads
        phix_stats = DECONTAMBAM_PHIX.out.stats
    }
    else {
        decontaminated_reads = reads
        phix_stats = Channel.empty()
    }

    decontaminated_reads = decontaminated_reads.map { meta, reads ->
        [meta + ['decontam_phix_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    // decontaminated_reads.view{ meta, _reads -> "decontaminated_phix_reads - [${meta.id}, ${meta.platform}, ${meta.single_end}] - ${meta.decontam_phix_read_count}" }
    decontaminated_reads = decontaminated_reads.filter { meta, _reads ->
        meta.decontam_phix_read_count > 0
    }

    if (host_genome != null) {

        // chunked_decontaminated_reads = decontaminated_reads
        //     .flatMap{ meta, fastas ->
        //         def chunks = []
        //         if(meta.single_end) {
        //             chunks = fastas.splitFasta(file: true, compress: true, size: 500.KB)
        //         } else {
        //             chunks = [
        //                 fastas[0].splitFasta(file: true, compress: true, size: 500.KB),
        //                 fastas[1].splitFasta(file: true, compress: true, size: 500.KB)
        //             ].transpose()
        //         }
        //         chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
        //     }
        // chunked_decontaminated_reads.view{ "chunked_decontaminated_reads - ${it}" }

        // SEQKIT_SPLIT2_HOST(decontaminated_reads)
        CHUNKFASTX_HOST(decontaminated_reads, file("${projectDir}/bin/chunk_fastx.py"))
        chunked_decontaminated_reads = CHUNKFASTX_HOST.out.reads
            .flatMap{ meta, chunks ->
                if (chunks instanceof List) {
                    def grouped_chunks = [:]
                    chunks.collect{
                        chunk ->
                        def part = (chunk.name =~ /\.(chunk-\d+)\./)[0][1]
                        if (grouped_chunks.containsKey(part)) {
                            grouped_chunks[part].add(chunk)
                        } else {
                            grouped_chunks[part] = [chunk]
                        }
                    }
                    return grouped_chunks.collect{ _part, chunk -> tuple(groupKey(meta, grouped_chunks.size()), tuple(*chunk)) }
                } else {
                    return [tuple(groupKey(meta, 1), chunks)]
                }
            }
        // chunked_decontaminated_reads.view{ "chunked_decontaminated_reads - ${it}" }

        BWAMEM2_ALIGN_HOST(
            chunked_decontaminated_reads,
            host_genome_index,
            host_genome_fasta,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_HOST.out.versions)

        COMBINEBAM_HOST(
            BWAMEM2_ALIGN_HOST.out.bam.groupTuple()
        )

        DECONTAMBAM_HOST(
            COMBINEBAM_HOST.out.concatenated_result.map{ meta, bam ->
                [meta, bam, meta.single_end==false, "short_read_host"]
            }
        )
        ch_versions = ch_versions.mix(DECONTAMBAM_HOST.out.versions)

        decontaminated_reads = DECONTAMBAM_HOST.out.unmapped_reads
        host_stats = DECONTAMBAM_HOST.out.stats
    }
    else {
        host_stats = Channel.empty()
    }

    decontaminated_reads = decontaminated_reads.map { meta, reads ->
        [meta + ['decontam_host_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    // decontaminated_reads.view{ meta, _reads -> "decontaminated_host_reads - [${meta.id}, ${meta.platform}, ${meta.single_end}] - ${meta.decontam_host_read_count}" }
    decontaminated_reads = decontaminated_reads.filter { meta, _reads ->
        meta.decontam_host_read_count > 0
    }

    emit:
    decontaminated_reads = decontaminated_reads
    host_stats           = host_stats
    phix_stats           = phix_stats
    versions             = ch_versions
}
