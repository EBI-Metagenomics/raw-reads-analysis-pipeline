/*
    ~~~~~~~~~~~~~~~~~~
     Imports
    ~~~~~~~~~~~~~~~~~~
*/
include { FETCHDB                 } from '../subworkflows/local/fetchdb/main'
include { QC                      } from '../subworkflows/local/qc/main'
include { READSMERGE              } from '../subworkflows/local/reads_qc/merge/main'
include { DECONTAM_SHORTREAD      } from '../subworkflows/local/decontam_shortread/main'
include { DECONTAM_LONGREAD       } from '../subworkflows/local/decontam_longread/main'

include { RRNA_EXTRACTION         } from '../subworkflows/ebi-metagenomics/rrna_extraction/main'
include { MAPSEQ_OTU_KRONA        } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { KRONA_KTIMPORTTEXT      } from '../modules/ebi-metagenomics/krona/ktimporttext/main'
include { MOTUS_PROFILE           } from '../modules/nf-core/motus/profile/main'
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                 } from '../modules/nf-core/multiqc/main'

include { samplesheetToList       } from 'plugin/nf-schema'


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
    db_ch = Channel.from(params.databases.collect{ k, v ->
        if((v instanceof Map) && (v.containsKey('base_dir'))){
            return [id: k] + v
        }
    }).filter{ it }
    // db_ch.view{ "db_ch - ${it}" }
    FETCHDB(db_ch, params.databases.cache_path)
    dbs_path_ch = FETCHDB.out.dbs
    // dbs_path_ch.view{ "dbs_path_ch - ${it}" }

    dbs_path_ch.branch{ meta, _fp ->
        motus: meta.id=='motus'
        host_genome: meta.id=='host_genome'
        host_genome_minimap2: meta.id=='host_genome_minimap2'
        phix: meta.id=='phix'
        rfam: meta.id=='rfam'
        silva_ssu: meta.id=='silva_ssu'
        silva_lsu: meta.id=='silva_lsu'
    }.set{ dbs }

    dbs.motus.view{ "dbs.motus - ${ it }" }
    dbs.host_genome.view{ "dbs.host_genome - ${ it }" }
    dbs.host_genome_minimap2.view{ "dbs.host_genome_minimap2 - ${ it }" }
    dbs.phix.view{ "dbs.phix - ${ it }" }
    dbs.rfam.view{ "dbs.rfam - ${ it }" }
    dbs.silva_ssu.view{ "dbs.silva_ssu - ${ it }" }
    dbs.silva_lsu.view{ "dbs.silva_lsu - ${ it }" }

    classified_reads = fetch_reads_transformed.map { meta, reads ->
        // Long reads
        if ( ["ont", "pb"].contains( meta.platform ) ) {
            return [meta + [long_reads: true], reads]
        // Short reads
        } else {
            return [meta + [short_reads: true], reads]
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

    DECONTAM_SHORTREAD(
        reads_to_analyse.short_reads,
        dbs.host_genome,
        dbs.phix
    )
    ch_versions = ch_versions.mix(DECONTAM_SHORTREAD.out.versions)

    // DECONTAM_LONGREAD(
    //     reads_to_analyse.long_reads,
    //     dbs.host_genome_minimap2,
    // )
    // ch_versions = ch_versions.mix(DECONTAM_LONGREAD.out.versions)

    // ch_qc_out = DECONTAM_SHORTREAD.out.mix(DECONTAM_LONGREAD.out)

    // FASTQC(
    //     ch_qc_out.qc_reads
    // )
    // ch_qc_out.view { "ch_qc_out - ${it}" }

    // // mOTUs
    // MOTUS_PROFILE(
    //     ch_qc_out,
    //     dbs['motus'],
    // )
    // ch_versions = ch_versions.mix(MOTUS_PROFILE.out.versions)
    // MOTUS_PROFILE.out.out.view { "motus.out - ${it}" }

    // // KRONA
    // // KRONA_KTIMPORTTEXT(MOTUS_PROFILE.out.out)
    // // ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    // // rrna_extraction
    // READSMERGE(ch_qc_out)
    // ch_versions = ch_versions.mix(READSMERGE.out.versions)

    // (db_meta, db_fp) = dbs['rfam']
    // RRNA_EXTRACTION(
    //     READSMERGE.out.reads_fasta,
    //     Channel.value("${db_fp}/${db_meta.files.ribosomal_models_file}"),
    //     Channel.value("${db_fp}/${db_meta.files.ribosomal_claninfo_file}")
    // )
    // ch_versions = ch_versions.mix(RRNA_EXTRACTION.out.versions)

    // (db_meta, db_fp) = dbs['silva_lsu']
    // lsu_db = [db_meta, ["${db_fp}/${db_meta.files.fasta}",
    //                     "${db_fp}/${db_meta.files.tax}",
    //                     "${db_fp}/${db_meta.files.otu}",
    //                     "${db_fp}/${db_meta.files.mscluster}"]

    // (db_meta, db_fp) = dbs['silva_ssu']
    // ssu_db = [db_meta, ["${db_fp}/${db_meta.files.fasta}",
    //                     "${db_fp}/${db_meta.files.tax}",
    //                     "${db_fp}/${db_meta.files.otu}",
    //                     "${db_fp}/${db_meta.files.mscluster}"]
    // lsu_ch = RRNA_EXTRACTION.out.lsu_fasta.combine(Channel.value(lsu_db))
    // ssu_ch = RRNA_EXTRACTION.out.ssu_fasta.combine(Channel.value(ssu_db))
    // rrna_ch = lsu_ch.mix(ssu_ch)
    // rrna_chs = rrna_ch.multiMap { seqs, db ->
    //     seqs: seqs
    //     db: db
    // }
    // MAPSEQ_OTU_KRONA(rrna_chs.seqs, rrna_chs.db)
    // ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)
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
