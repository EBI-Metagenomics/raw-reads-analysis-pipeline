include { INFERNAL_CMSEARCH           } from '../../../modules/ebi-metagenomics/infernal/cmsearch/main'
include { CMSEARCHTBLOUTDEOVERLAP     } from '../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main'
include { EASEL_ESLSFETCH             } from '../../../modules/ebi-metagenomics/easel/eslsfetch/main'
include { EXTRACTCOORDS               } from '../../../modules/ebi-metagenomics/extractcoords/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'

workflow RRNA_EXTRACTION {

    take:
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    rfam         // file: rfam for cmsearch
    claninfo     // file: claninfo for cmsearchtbloutdeoverlap

    main:
    ch_versions = Channel.empty()

    ch_fasta = ch_fasta.map { meta, reads ->
        [meta + ['fasta_read_count': reads.countFasta()], reads]
    }
    // ch_fasta.view{ meta, _reads -> "ch_fasta - [${meta.id}, ${meta.platform}, ${meta.single_end}] - ${meta.fasta_read_count}" }
    // decontaminated_reads = decontaminated_reads.filter { meta, _reads ->
    //     meta.decontam_host_read_count > 0
    // }

    ch_chunked_fasta = ch_fasta
        .flatMap{ meta, fasta ->
            def chunks = fasta.splitFasta(file: true, size: 1.MB)
            chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
        }

    INFERNAL_CMSEARCH(
        ch_chunked_fasta,
        rfam
    )
    ch_versions = ch_versions.mix(INFERNAL_CMSEARCH.out.versions.first())

    COMBINEHMMSEARCHTBL(
        INFERNAL_CMSEARCH.out.cmsearch_tbl.groupTuple()
    )

    CMSEARCHTBLOUTDEOVERLAP(
        COMBINEHMMSEARCHTBL.out.concatenated_result,
        claninfo
    )
    ch_versions = ch_versions.mix(CMSEARCHTBLOUTDEOVERLAP.out.versions.first())

    ch_easel = ch_fasta
                .join(CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped)

    EASEL_ESLSFETCH(ch_easel)
    ch_versions = ch_versions.mix(EASEL_ESLSFETCH.out.versions.first())

    EXTRACTCOORDS(
        EASEL_ESLSFETCH.out.easel_coords,
        EASEL_ESLSFETCH.out.matched_seqs_with_coords
    )
    ch_versions = ch_versions.mix(EXTRACTCOORDS.out.versions.first())

    emit:
    cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped   // channel: [ val(meta), [ deoverlapped ] ]
    easel_out = EASEL_ESLSFETCH.out.easel_coords                                        // channel: [ val(meta), [ fasta ] ]
    ssu_fasta = EXTRACTCOORDS.out.ssu_fasta                                             // channel: [ val(meta), [ fasta ] ]
    lsu_fasta = EXTRACTCOORDS.out.lsu_fasta                                             // channel: [ val(meta), [ fasta ] ]
    concat_ssu_lsu_coords = EXTRACTCOORDS.out.concat_ssu_lsu_coords                     // channel: [ val(meta), [ txt ] ]
    versions = ch_versions                                                              // channel: [ versions.yml ]
}

