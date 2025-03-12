/*
    ~~~~~~~~~~~~~~~~~~
     Imports
    ~~~~~~~~~~~~~~~~~~
*/
// include { FETCHREADS              } from '../subworkflows/local/fetchreads/main'
// include { FETCHDB                 } from '../subworkflows/local/fetchdb/main'
include { QC                      } from '../subworkflows/local/qc/main'
include { READSMERGE              } from '../subworkflows/local/reads_qc/merge/main'
include { DECONTAM_SHORTREAD      } from '../subworkflows/local/decontam_shortread/main'
include { DECONTAM_LONGREAD       } from '../subworkflows/local/decontam_longread/main'

include { KRONA_KTIMPORTTEXT      } from '../modules/ebi-metagenomics/krona/ktimporttext/main'
include { MOTUS_PROFILE           } from '../modules/nf-core/motus/profile/main'
include { RRNA_EXTRACTION         } from '../subworkflows/ebi-metagenomics/rrna_extraction/main'
include { MAPSEQ_OTU_KRONA        } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
// include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// Import samplesheetToList from nf-schema //
include { samplesheetToList       } from 'plugin/nf-schema'

include { FASTQC } from '../modules/nf-core/fastqc/main'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow PIPELINE {
    main:
    ch_versions = Channel.empty()

    def groupReads = { study_accession, reads_accession, fq1, fq2, library_layout, library_strategy, platform -> [
        ['id': reads_accession,
         'study_accession': study_accession,
         'single_end': (fq2 == []) ? true : false,
         'library_layout': library_layout,
         'library_strategy': library_strategy,
         'platform': params.platform ?: platform,],
        (fq2 == []) ? [fq1] : [fq1, fq2],
    ]}
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))

    // [ study, sample, read1, [read2], library_layout, library_strategy, platform]
    fetch_reads_transformed = samplesheet.map(groupReads)

    // Fetch databases
    dbs = params.databases.collectEntries { db_id, db_data ->
        if ((db_data instanceof Map) && db_data.containsKey('name')) {
            return [db_id, [meta: [id: db_id] + db_data,
                            file: file(db_data.path, checkIfExists: true)]]
        }
    }

    classified_reads = fetch_reads_transformed.map { meta, reads ->
        // Long reads //
        if ( ["ont", "pb"].contains( meta.platform ) ) {
            return [ meta + [long_reads: true], reads]
        // Short reads //
        } else {
            return [ meta + [short_reads: true], reads]
        }
    }

    // QC
    QC(classified_reads)
    ch_versions = ch_versions.mix(QC.out.versions)

    // DECONTAMINATION
    QC.out.qc_reads.branch { meta, _reads ->
            short_reads: meta.short_reads
            long_reads: meta.long_reads
        }
    .set { reads_to_analyse }

    db = dbs['host_genome']
    host_genome_db = [db.meta, "${db.file}/${db_meta.files.bwa_index_prefix}"]
    (db_meta, db_fp) = dbs['human_phix']
    human_phix_db = [db_meta, "${db_fp}/${db_meta.files.bwa_index_prefix}"]

    DECONTAM_SHORTREAD(
        reads_to_analyse.short_reads,
        Channel.value(host_genome_db),
        Channel.value(human_phix_db)
    )
    ch_versions = ch_versions.mix(DECONTAM_SHORTREAD.out.versions)

    db_meta, db_fp = dbs['host_genome_minimap2']
    // host_genome_minimap2_db = [db_meta, "${db_fp}/${db_meta.files.index}"]
    host_genome_minimap2_db = [db_meta, db_fp]
    DECONTAM_LONGREAD(
        reads_to_analyse.long_reads,
        Channel.value(host_genome_minimap2_db),
    )
    ch_versions = ch_versions.mix(LONG_READ_QC.out.versions)

    ch_qc_out = SHORT_READS_QC.out.mix(LONG_READ_QC.out)

    FASTQC(
        ch_qc_out.qc_reads
    )
    ch_qc_out.view { "ch_qc_out - ${it}" }

    // mOTUs
    MOTUS_PROFILE(
        ch_qc_out,
        dbs['motus'],
    )
    ch_versions = ch_versions.mix(MOTUS_PROFILE.out.versions)
    MOTUS_PROFILE.out.out.view { "motus.out - ${it}" }

    // KRONA
    // KRONA_KTIMPORTTEXT(MOTUS_PROFILE.out.out)
    // ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    // rrna_extraction
    READSMERGE(ch_qc_out)
    ch_versions = ch_versions.mix(READSMERGE.out.versions)

    (db_meta, db_fp) = dbs['rfam']
    RRNA_EXTRACTION(
        READSMERGE.out.reads_fasta,
        Channel.value("${db_fp}/${db_meta.files.ribosomal_models_file}"),
        Channel.value("${db_fp}/${db_meta.files.ribosomal_claninfo_file}")
    )
    ch_versions = ch_versions.mix(RRNA_EXTRACTION.out.versions)

    (db_meta, db_fp) = dbs['silva_lsu']
    lsu_db = [db_meta, ["${db_fp}/${db_meta.files.fasta}",
                        "${db_fp}/${db_meta.files.tax}",
                        "${db_fp}/${db_meta.files.otu}",
                        "${db_fp}/${db_meta.files.mscluster}"]

    (db_meta, db_fp) = dbs['silva_ssu']
    ssu_db = [db_meta, ["${db_fp}/${db_meta.files.fasta}",
                        "${db_fp}/${db_meta.files.tax}",
                        "${db_fp}/${db_meta.files.otu}",
                        "${db_fp}/${db_meta.files.mscluster}"]
    lsu_ch = RRNA_EXTRACTION.out.lsu_fasta.combine(Channel.value(lsu_db))
    ssu_ch = RRNA_EXTRACTION.out.ssu_fasta.combine(Channel.value(ssu_db))
    rrna_ch = lsu_ch.mix(ssu_ch)
    rrna_chs = rrna_ch.multiMap { seqs, db ->
        seqs: seqs
        db: db
    }
    MAPSEQ_OTU_KRONA(rrna_chs.seqs, rrna_chs.db)
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)
    // MULTIQC(READSTRIM.out.fastp_json.join(MOTUS.out.map{ [it[0],it[3]] }))

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'ASSEMBLY_ANALYSIS_PIPELINE_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { collated_versions }

    emit:
    versions = ch_versions
    collated_versions = collated_versions
}
