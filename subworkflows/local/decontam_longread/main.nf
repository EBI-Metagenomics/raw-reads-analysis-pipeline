include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { DECONTAMBAM    } from '../../../modules/local/decontambam/main'
include { COMBINEBAM } from '../../../modules/local/combinebam/main'
// include { SEQKIT_SPLIT2 } from '../../../modules/nf-core/seqkit/split2/main'
include { CHUNKFASTX } from '../../../modules/local/chunkfastx/main'
include { GZIPALL } from '../../../modules/local/gzipall/main'

workflow DECONTAM_LONGREAD {
    take:
    input_reads       // [ val(meta), path(reads) ]
    reference_genome  // [ val(meta2), path(reference_genome) ]

    main:
    def ch_versions = Channel.empty()

    reference_genome_index = reference_genome
        .map{ meta, fp -> [
            meta,
            file("${fp}/${meta.base_dir}/${meta.files.index}")
        ] }
        .first()

    // chunked_reads = input_reads
    //     .flatMap{ meta, fastas ->
    //         def chunks = fastas.splitFasta(file: true, compress: true, size: 500.KB)
    //         chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
    //     }
    // chunked_reads.view{ "chunked_reads - ${it}" }

    // SEQKIT_SPLIT2(input_reads)
    CHUNKFASTX(
        input_reads,
        file("${projectDir}/bin/chunk_fastx.py")
    )
    chunked_reads = CHUNKFASTX.out.reads
        .flatMap{ meta, chunks ->
            if (chunks instanceof List) {
                return chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
            } else {
                return [tuple(groupKey(meta, 1), chunks)]
            }
        }
    chunked_reads = chunked_reads.map{ meta, reads -> 
                                       [meta, reads, reads[0].name.endsWith('.gz')] } 
        .branch{ meta, reads, zip -> 
            to_zip: !zip
                return [meta, reads]
            already_zip: zip
                return [meta, reads] 
        } 
    GZIPALL(chunked_reads.to_zip)
    chunked_reads = chunked_reads.already_zip.mix(GZIPALL.out.files)
    // chunked_reads.view{ "chunked_reads - ${it}" }

    MINIMAP2_ALIGN(
        chunked_reads,
        reference_genome_index,
        true,
        false,
        "bai",
        false,
        false,
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    COMBINEBAM(
        MINIMAP2_ALIGN.out.bam.groupTuple()
    )

    DECONTAMBAM(
        COMBINEBAM.out.concatenated_result.map{ meta, bam ->
            [meta, bam, meta.single_end==false, "long_read_host"]
        }
    )
    ch_versions = ch_versions.mix(DECONTAMBAM.out.versions)

    decontaminated_reads = DECONTAMBAM.out.unmapped_reads.map { meta, reads ->
        [meta + ['decontam_host_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    // decontaminated_reads.view{ meta, _reads -> "decontaminated_host_reads - [${meta.id}, ${meta.platform}, ${meta.single_end}] - ${meta.decontam_host_read_count}" }
    decontaminated_reads = decontaminated_reads.filter { meta, _reads ->
        meta.decontam_host_read_count > 0
    }

    emit:
    decontaminated_reads = decontaminated_reads
    stats                = DECONTAMBAM.out.stats
    versions             = ch_versions
}
