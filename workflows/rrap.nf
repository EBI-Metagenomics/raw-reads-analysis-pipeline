/*
    ~~~~~~~~~~~~~~~~~~
     Imports
    ~~~~~~~~~~~~~~~~~~
*/
include { ACCESSION2SAMPLESHEET  } from '../modules/local/accession2samplesheet/main'  // JOINSAMPLESHEETS 
include { READSDECONTAM          } from "../modules/local/readsdecontam/main"
include { FETCHREADS             } from '../subworkflows/local/fetchreads/main'
include { FETCHDB                } from "../subworkflows/local/fetchdb/main"
include { READSQC                } from '../subworkflows/local/reads_qc/qc/main'
include { READSTRIM              } from '../subworkflows/local/reads_qc/trim/main'
include { READSMERGE             } from '../subworkflows/local/reads_qc/merge/main'

include { KRONA_KTIMPORTTEXT     } from '../modules/ebi-metagenomics/krona/ktimporttext/main'
include { MOTUS_PROFILE          } from '../modules/nf-core/motus/profile/main'
include { RRNA_EXTRACTION        } from "../subworkflows/ebi-metagenomics/rrna_extraction/main"
include { MAPSEQ_OTU_KRONA       } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
// include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// Import samplesheetToList from nf-schema //
include { samplesheetToList      } from 'plugin/nf-schema'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE {
    versions = Channel.empty()

    // Read input samplesheet and validate it using schema_input.json //
    if((params.input instanceof String) && (params.input.size()>0)) {
        samplesheet_fp = Channel.value(params.input)
    }else{
        samplesheet_fp = ACCESSION2SAMPLESHEET(params.accession)
    }
   
    samplesheet_multi = samplesheet_fp
        .map{ samplesheetToList(it, "${workflow.projectDir}/assets/schema_input.json") }
        .multiMap{
            meta: it.collect{ it[0] }
            fp1: it.collect{ it[1] }
            fp2: it.collect{ it[2] }
        }
    
    samplesheet = samplesheet_multi.meta.flatten()
           .merge(samplesheet_multi.fp1.flatten())
           .merge(samplesheet_multi.fp2.flatten())
    // samplesheet.view{ "samplesheet - ${it}" }

    // Fetch reads from samplesheet
    fetch_ch_split = samplesheet.multiMap{ meta, fp1, fp2 ->
        fp1: [[id: meta.id, accession: meta.accession, read: 1], fp1]
        fp2: [[id: meta.id, accession: meta.accession, read: 2], fp2] 
    }
    fetch_ch = fetch_ch_split.fp1.mix( fetch_ch_split.fp2 )
        .toSortedList{ a, b -> "${a[0]['id']}-${a[0]['read']}" <=> "${b[0]['id']}-${b[0]['read']}" }
        .flatMap()
        .filter{ meta, fp -> fp ? true: false }
        .map{ meta, fp -> [meta.id, [meta, fp]] }

    sizes = fetch_ch.groupTuple().map{ meta_id, reads_rows -> 
        [meta_id, reads_rows.collect{ meta, reads -> if(reads){ return reads } }.size()] 
    }

    fetch_ch = fetch_ch.combine(sizes, by: 0).map{ 
        meta_id, reads_row, size -> 
        def (meta, reads) = [reads_row[0], reads_row[1]]
        [meta + [size: size], reads] 
    }
    // fetch_ch.view{ "fetch_ch - ${it}" }
    FETCHREADS(fetch_ch, Channel.value(params.reads_cache_path))
    // FETCHREADS.out.view{ "FETCHREADS.out - ${it}" }
     
    // Fetch databases 
    db_ch = Channel.from(params.databases.collect{ k,v -> 
        if((v instanceof Map) && (v.containsKey('name'))){ 
            return [[id: k], v] 
        } 
    }).filter{ it }
    // db_ch.view{ "db_ch - ${it}" } 
    FETCHDB(db_ch, Channel.value(params.databases.cache_path))
    dbs = FETCHDB.out.dbs
    // dbs.view{ "dbs - ${it}" }
  
    //QC  
    qc_ch = FETCHREADS.out
        .map{ meta,fp -> 
              def meta_ = [id: meta.id, accession: meta.accession, size: meta.size]
              [groupKey(meta_,meta_.size), fp] }
        .groupTuple() 
        .map{ meta,fps -> 
              [[id: meta.id, 
                accession: meta.accession, 
                single_end: meta.size>1 ? false:true], 
               fps]}
    // qc_ch.view{ "qc_ch - ${it}" }
    
    READSQC(qc_ch, false)
    versions = versions.mix(READSQC.out.versions)
    READSTRIM(READSQC.out.passed_reads)
    versions = versions.mix(READSTRIM.out.versions)
    // READSDECONTAM
    READSDECONTAM(
        READSTRIM.out.reads, 
        dbs.filter{ meta, db -> meta.id == 'host_genome' }.first()
    )
    // versions = versions.mix(READSDECONTAM.out.versions)
    READSDECONTAM.out.view{ "READSDECONTAM.out - ${it}" }    
 
    // mOTUs
    MOTUS_PROFILE(
        READSDECONTAM.out, 
        dbs.filter{ meta, db -> meta.id == 'motus' }.map{ meta,db -> db }.first()
    )
    versions = versions.mix(MOTUS_PROFILE.out.versions)
    MOTUS_PROFILE.out.out.view{ "motus.out - ${it}" }
   
    // KRONA
    // KRONA_KTIMPORTTEXT(MOTUS_PROFILE.out.out)
    // versions = versions.mix(KRONA_KTIMPORTTEXT.out.versions) 
 
    // rrna_extraction
    READSMERGE(READSDECONTAM.out)
    versions = versions.mix(READSMERGE.out.versions)

    rfam_db = dbs.filter{ meta, db -> meta.id == 'rfam' }
    RRNA_EXTRACTION(
        READSMERGE.out.reads_fasta, 
        rfam_db.map{ meta, db -> "${db}/${params.databases.rfam.files.ribosomal_models_file}" }.first(), 
        rfam_db.map{ meta, db -> "${db}/${params.databases.rfam.files.ribosomal_claninfo_file}" }.first()
    )    
    versions = versions.mix(RRNA_EXTRACTION.out.versions)
    rrna_dbs = rfam_db
        .filter{ meta,db -> meta.id in ['silva_ssu','silva_lsu']} 
        .map{ meta,db -> [meta, [
            "${db}/${params.databases[meta.id].files.fasta}",
            "${db}/${params.databases[meta.id].files.tax}",
            "${db}/${params.databases[meta.id].files.otu}",
            "${db}/${params.databases[meta.id].files.mscluster}",
            meta.id,
        ]] }
    lsu_db = rrna_dbs.filter{ meta, db -> meta.id == 'silva_lsu' }
                     .map{ meta, db -> db } 
                     .first()
    ssu_db = rrna_dbs.filter{ meta, db -> meta.id == 'silva_ssu' }
                     .map{ meta, db -> db }
                     .first() 
    lsu_ch = RRNA_EXTRACTION.out.lsu_fasta.combine(lsu_db)
    ssu_ch = RRNA_EXTRACTION.out.ssu_fasta.combine(ssu_db)
    rrna_ch = lsu_ch.mix(ssu_ch)
    rrna_chs = rrna_ch.multiMap{ seqs, db ->
        seqs: seqs
        db: db
    }
    MAPSEQ_OTU_KRONA(rrna_chs.seqs, rrna_chs.db) 
    versions = versions.mix(MAPSEQ_OTU_KRONA.out.versions)
    // MULTIQC(READSTRIM.out.fastp_json.join(MOTUS.out.map{ [it[0],it[3]] }))
   
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'ASSEMBLY_ANALYSIS_PIPELINE_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { collated_versions }
 
    emit:
    versions = versions
}
