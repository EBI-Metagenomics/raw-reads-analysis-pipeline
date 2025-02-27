/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { QC } from '../subworkflows/qc_swf'
include { MAPSEQ_OTU_KRONA } from '../subworkflows/mapseq_otu_krona_swf'
include { CMSEARCH_SUBWF } from '../subworkflows/cmsearch_swf'
include { FETCHTOOL_RAWREADS } from '../modules/fetchtool'
include { MOTUS } from '../modules/motus'
include { MULTIQC } from '../modules/multiqc'
include { FETCH_READS } from '../subworkflows/fetch_reads'

/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_MOTUS_DB } from '../subworkflows/prepare_dbs'
include { DOWNLOAD_HOST_REFERENCE_GENOME } from '../subworkflows/prepare_dbs'
include { DOWNLOAD_RFAM } from '../subworkflows/prepare_dbs'
include { DOWNLOAD_MAPSEQ_SSU } from '../subworkflows/prepare_dbs'
include { DOWNLOAD_MAPSEQ_LSU } from '../subworkflows/prepare_dbs'

// Import samplesheetToList from nf-schema //
include { samplesheetToList } from 'plugin/nf-schema'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE {
    // Read input samplesheet and validate it using schema_input.json //
    samplesheet = samplesheetToList(params.input, "./assets/schema_input.json") 
    // samplesheet_ch = Channel.fromList(samplesheet)
    // samplesheet_ch.view{ "samplesheet - ${it}" }

    sample2fp_list = []
    samplesheet.each{
        meta, fq1, fq2, fqb, fq1_md5, fq2_md5, fqb_md5 -> 
        if(fq1) {
            sample2fp_list.add([meta,fq1,fq1_md5,'reads'])
        }
        if(fq2) {
            sample2fp_list.add([meta,fq2,fq2_md5,'reads'])
        }
        if(fqb) {
            sample2fp_list.add([meta,fqb,fqb_md5,'barcodes'])
        }
    }
    sizes = [:]
    sample2fp_list.each{
        meta, fp, md5, t ->
        if(sizes[meta.id]) {
            sizes[meta.id] += 1
        }else{
            sizes[meta.id] = 1
        }
    }

    fetch_ch = Channel.fromList(sample2fp_list)
    // fetch_ch.view{ "fetch_ch - ${it}"}
    FETCH_READS(fetch_ch)
    // FETCH_READS.out.view{ "FETCH_READS.out - ${it}" }
     
    qc_ch = FETCH_READS.out.map{ 
        meta,fp,t -> 
        [groupKey(meta,sizes[meta.id]),[fp,t]] 
    }.groupTuple() 
    qc_ch = qc_ch.map{ 
        meta,fps -> 
        def fps_d = [:]
        fps.sort{ a,b -> a[0] <=> b[0] }.each{ fp,t -> 
            if(fps_d[t]) {
                fps_d[t].add(fp)
            }else{
                fps_d[t] = [fp]
            }
        }
        return [meta,fps_d]
    }
    qc_ch = qc_ch.filter{ 
        meta,fps -> 
        (fps.reads) && (fps.reads.size()>0)
    }.map{ 
        meta,fps -> 
        [meta, [reads: fps.reads, 
                mode: fps.reads.size()>1 ? 'paired':'single']]
    }
    // qc_ch.view{ "qc_ch - ${it}" }

    DOWNLOAD_HOST_REFERENCE_GENOME()
    host_genome_db = DOWNLOAD_HOST_REFERENCE_GENOME.out.host_genome_db
    // host_genome_db.view{ "host_genome_db - ${it}"}

    QC(
        qc_ch,
        host_genome_db,
    )
    QC.out.merged_reads.view{ "QC.out.merged_reads - ${it}" }


    // mOTUs
    DOWNLOAD_MOTUS_DB()
    motus_db = DOWNLOAD_MOTUS_DB.out
    motus_db.view{ "motus_db - ${it}" }

    MOTUS(QC.out.merged_reads, motus_db)
    MOTUS.out.view{ "MOTUS.out - ${it}" }
    
    
    // cmsearch 
    DOWNLOAD_RFAM()
    rfam_dbs = DOWNLOAD_RFAM.out.rfam_db_dir

    CMSEARCH_SUBWF(QC.out.sequence, rfam_dbs)
    
     
    // mapseq with silva ssu and lsu
    DOWNLOAD_MAPSEQ_LSU()
    mapseq_lsu_db_dir = DOWNLOAD_MAPSEQ_LSU.out.mapseq_lsu_db_dir
    mapseq_lsu_db_dir.view{ "mapseq_lsu_db_dir - ${it}" }
    
    DOWNLOAD_MAPSEQ_SSU()
    mapseq_ssu_db_dir = DOWNLOAD_MAPSEQ_SSU.out.mapseq_ssu_db_dir
    mapseq_ssu_db_dir.view{ "mapseq_ssu_db_dir - ${it}" }
    
    lsu_ch = CMSEARCH_SUBWF.out.cmsearch_lsu_fasta.combine(mapseq_lsu_db_dir)
    ssu_ch = CMSEARCH_SUBWF.out.cmsearch_ssu_fasta.combine(mapseq_ssu_db_dir)
    mapseq_in = lsu_ch.mix(ssu_ch)
    mapseq_in = mapseq_in.map{ meta_seq,seq,meta_db,db ->
        [[seq_id: meta_seq.id, db_id: meta_db.id], seq, db] }
    mapseq_in.view{ "mapseq_in - ${it}" }
    MAPSEQ_OTU_KRONA(mapseq_in)
    
    MULTIQC(QC.out.fastp_json.join(MOTUS.out.map{ [it[0],it[3]] }))
}
