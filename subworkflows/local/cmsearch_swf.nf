/*
 * Cmsearch and deoverlap-cmsearch
 */

include { CMSEARCH } from '../modules/cmsearch'
include { CMSEARCH_DEOVERLAP } from '../modules/cmsearch_deoverlap'
include { EASEL_EXTRACT_BY_COORD } from '../modules/easel'
include { EXTRACT_MODELS } from '../modules/extract_coords'

process COLLATE_FILES {

    publishDir(
        "${params.outdir}/cmsearch/",
        mode: 'copy'
    )

    container 'quay.io/biocontainers/infernal:1.1.4--pl5321hec16e2b_1'
    tag "${meta.id}"
    label 'process_single'

    stageInMode 'copy'

    input:
    tuple val(meta), file(cmsearch), file(deoverlap)

    output:
    tuple val(meta), file("${meta.id}_matched_seqs_with_coords.tbl"), file("${meta.id}_matched_seqs_with_coords_deoverlap.tbl")

    script:
    """
    head -n 1 ${cmsearch} > "${meta.id}_matched_seqs_with_coords.tbl"
    grep -v '^#' ${cmsearch} | grep . >> "${meta.id}_matched_seqs_with_coords.tbl"

    head -n 1 ${cmsearch} > "${meta.id}_matched_seqs_with_coords_deoverlap.tbl"
    cat ${deoverlap} >> "${meta.id}_matched_seqs_with_coords_deoverlap.tbl"
    """
}

workflow CMSEARCH_SUBWF {
    take:
        sequences
        rfam_dbs
    main:
        rfam_dbs_split = rfam_dbs.multiMap{ meta,d ->
            def db_fps = params.databases.rfam.files
            ribo_models: [meta, file("${d}/${db_fps.ribosomal_models_file}")]
            other_models: [meta, file("${d}/${db_fps.other_models_file}")]
            ribo_claninfo: [meta, file("${d}/${db_fps.ribosomal_claninfo_file}")]
            other_claninfo: [meta, file("${d}/${db_fps.other_claninfo_file}")]
        } 

        // run
        CMSEARCH(sequences, rfam_dbs_split.ribo_models)
        CMSEARCH_DEOVERLAP(CMSEARCH.out, rfam_dbs_split.ribo_claninfo) 
        
        // save outputs
        CMSEARCH.out.collectFile(name: "cmsearch.tbl", newLine: true)
        CMSEARCH_DEOVERLAP.out.collectFile(name: "deoverlapped.tbl")

        EASEL_EXTRACT_BY_COORD(sequences.join(CMSEARCH_DEOVERLAP.out))
        EXTRACT_MODELS(EASEL_EXTRACT_BY_COORD.out)
        
        extract_models_out = EXTRACT_MODELS.out.multiMap{ meta,a,b,c ->
            directory: [meta, a]
            ssu_fasta: [meta, b]
            lsu_fasta: [meta, c]
        }
        
        output_ch = CMSEARCH.out.join(CMSEARCH_DEOVERLAP.out) 
        COLLATE_FILES(output_ch)
    emit:
        cmsearch_lsu_fasta = extract_models_out.lsu_fasta
        cmsearch_ssu_fasta = extract_models_out.ssu_fasta
        seq_cat = extract_models_out.directory
}
