/*
    ~~~~~~~~~~~~~~~~~~
     Imports
    ~~~~~~~~~~~~~~~~~~
*/
include { FETCHDB } from '../subworkflows/local/fetchdb/main'
include { QC } from '../subworkflows/local/qc/main'
include { READSMERGE } from '../subworkflows/local/readsmerge/main'
include { DECONTAM_SHORTREAD } from '../subworkflows/local/decontam_shortread/main'
include { DECONTAM_LONGREAD } from '../subworkflows/local/decontam_longread/main'
include { MOTUS_KRONA } from '../subworkflows/local/motus_krona/main'
include { SEQTK_SEQ } from '../modules/ebi-metagenomics/seqtk/seq/main'

include { RRNA_EXTRACTION } from '../subworkflows/local/rrna_extraction/main'
include { MAPSEQ_OTU_KRONA } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { PROFILE_HMMSEARCH_PFAM } from '../subworkflows/local/profile_hmmsearch_pfam/main'
include { samplesheetToList } from 'plugin/nf-schema'


workflow PIPELINE {
    main:
    ch_versions = Channel.empty()

    // Fetch databases
    db_ch = Channel
        .from(
            params.databases.collect { k, v ->
                if ((v instanceof Map) && v.containsKey('base_dir')) {
                    return [id: k] + v
                }
            }
        )
        .filter { it }
    // db_ch.view{ "db_ch - ${it}" }
    FETCHDB(db_ch, params.databases.cache_path)
    dbs_path_ch = FETCHDB.out.dbs
    // dbs_path_ch.view{ "dbs_path_ch - ${it}" }

    if (!params.download_dbs){

    dbs_path_ch
        .branch { meta, _fp ->
            motus: meta.id == 'motus'
            host_genome: meta.id == 'host_genome'
            host_genome_minimap2: meta.id == 'host_genome_minimap2'
            phix: meta.id == 'phix'
            rfam: meta.id == 'rfam'
            silva_ssu: meta.id == 'silva_ssu'
            silva_lsu: meta.id == 'silva_lsu'
            pfam: meta.id == 'pfam'
        }
        .set { dbs }

    // dbs.motus.view{ "dbs.motus - ${ it }" }
    // dbs.host_genome.view{ "dbs.host_genome - ${ it }" }
    // dbs.host_genome_minimap2.view{ "dbs.host_genome_minimap2 - ${ it }" }
    // dbs.phix.view{ "dbs.phix - ${ it }" }
    // dbs.rfam.view{ "dbs.rfam - ${ it }" }
    // dbs.silva_ssu.view{ "dbs.silva_ssu - ${ it }" }
    // dbs.silva_lsu.view{ "dbs.silva_lsu - ${ it }" }
    // dbs.pfam.view{ "dbs.pfam - ${ it }" }

    // Parse samplesheet and fetch reads
    def groupReads = { study_accession, reads_accession, fq1, fq2, library_layout, library_strategy, platform ->
        [
            [
                'id': reads_accession,
                'study_accession': study_accession,
                'single_end': fq2 == [] ? true : false,
                'library_layout': library_layout,
                'library_strategy': library_strategy,
                'platform': platform
            ],
            fq2 == [] ? file(fq1) : [file(fq1), file(fq2)]
        ]
    }
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))

    // [ study, sample, read1, [read2], library_layout, library_strategy, platform]
    fetch_reads_transformed = samplesheet.map(groupReads)

    classified_reads = fetch_reads_transformed.map { meta, reads ->
        // Long reads
        if (["ont", "pb"].contains(meta.platform)) {
            return [meta + [long_reads: true], reads]
        }
        else {
            return [meta + [short_reads: true], reads]
        }
    }

    // Get read count per fastq row
    classified_reads = classified_reads.map { meta, reads ->
        [meta + ['read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    // classified_reads.view{ meta, _reads -> "classified_reads - ${meta.id} - ${meta.read_count}" }
    classified_reads = classified_reads.filter { meta, _reads ->
        meta.read_count > 0
    }

    // QC
    if (params.skip_qc) {
        classified_reads.set{ qc_reads }
	qc_stats = Channel.empty()
    }
    else {
        QC(classified_reads)
        ch_versions = ch_versions.mix(QC.out.versions)

        qc_reads = QC.out.fastq
        qc_stats = QC.out.fastp_json
    }


    // Get read count per fastq row
    qc_reads = qc_reads.map { meta, reads ->
        [meta + ['qc_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    // qc_reads.view{ meta, _reads -> "qc_reads - ${meta.id} - ${meta.qc_read_count}" }
    qc_reads = qc_reads.filter { meta, _reads ->
        meta.qc_read_count > 0
    }

    // DECONTAMINATION
    if (params.skip_decontam) {
        qc_reads.set { clean_reads }
	decontam_stats = Channel.empty()
    }
    else {
        qc_reads
            .branch { meta, _reads ->
                short_reads: meta.short_reads
                long_reads: meta.long_reads
            }
            .set { reads_to_analyse }

        DECONTAM_SHORTREAD(
            reads_to_analyse.short_reads,
            dbs.host_genome,
            dbs.phix,
        )
        ch_versions = ch_versions.mix(DECONTAM_SHORTREAD.out.versions)

        DECONTAM_LONGREAD(
            reads_to_analyse.long_reads,
            dbs.host_genome_minimap2,
        )
        ch_versions = ch_versions.mix(DECONTAM_LONGREAD.out.versions)

        clean_reads = DECONTAM_SHORTREAD.out.decontaminated_reads
            .mix(DECONTAM_LONGREAD.out.decontaminated_reads)

        decontam_stats = DECONTAM_SHORTREAD.out.phix_stats
            .mix(DECONTAM_SHORTREAD.out.host_stats)
            .mix(DECONTAM_LONGREAD.out.stats)
    }

    // // Get read count per fastq row
    // clean_reads = clean_reads.map { meta, reads ->
    //     [meta + ['clean_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    // }
    // clean_reads.view{ meta, _reads -> "clean_reads - [${meta.id}, ${meta.platform}, ${meta.single_end}] - ${meta.clean_read_count}" }
    // clean_reads = clean_reads.filter { meta, _reads ->
    //     meta.clean_read_count > 0
    // }

    // FASTQC(clean_reads)
    // ch_versions = ch_versions.mix(FASTQC.out.versions)

    motus_db = dbs.motus
        .map { meta, fp ->
            file("${fp}/${meta.base_dir}")
        }
        .first()
    // motus_db.view{ "motus_db - ${it}" }

    // mOTUs
    MOTUS_KRONA(
        clean_reads.map { meta, reads ->
            [meta, meta.single_end ? [reads] : reads]
        },
        motus_db
    )
    ch_versions = ch_versions.mix(MOTUS_KRONA.out.versions)

    // rrna_extraction
    READSMERGE(clean_reads)
    ch_versions = ch_versions.mix(READSMERGE.out.versions)

    rfam_db = dbs.rfam
        .map { meta, fp ->
            file("${fp}/${meta.base_dir}/${meta.files.ribosomal_models_file}")
        }
        .first()
    // rfam_db.view{ "rfam_db - ${it}" }

    claninfo_db = dbs.rfam
        .map { meta, fp ->
            file("${fp}/${meta.base_dir}/${meta.files.ribosomal_claninfo_file}")
        }
        .first()
    // claninfo_db.view{ "claninfo_db - ${it}" }

    RRNA_EXTRACTION(
        READSMERGE.out.reads_fasta,
        rfam_db,
        claninfo_db,
    )
    ch_versions = ch_versions.mix(RRNA_EXTRACTION.out.versions)

    lsu_db = dbs.silva_lsu
        .map { meta, fp ->
            [
                [
                    file("${fp}/${meta.base_dir}/${meta.files.fasta}"),
                    file("${fp}/${meta.base_dir}/${meta.files.tax}"),
                    file("${fp}/${meta.base_dir}/${meta.files.otu}"),
                    file("${fp}/${meta.base_dir}/${meta.files.mscluster}"),
                    meta.id,
                ]
            ]
        }
        .first()
    // lsu_db.view{ "lsu_db - ${it}" }

    ssu_db = dbs.silva_ssu
        .map { meta, fp ->
            [
                [
                    file("${fp}/${meta.base_dir}/${meta.files.fasta}"),
                    file("${fp}/${meta.base_dir}/${meta.files.tax}"),
                    file("${fp}/${meta.base_dir}/${meta.files.otu}"),
                    file("${fp}/${meta.base_dir}/${meta.files.mscluster}"),
                    meta.id,
                ]
            ]
        }
        .first()
    // ssu_db.view{ "ssu_db - ${it}" }

    lsu_ch = RRNA_EXTRACTION.out.lsu_fasta.combine(lsu_db)
    ssu_ch = RRNA_EXTRACTION.out.ssu_fasta.combine(ssu_db)
    rrna_ch = lsu_ch.mix(ssu_ch)
    rrna_chs = rrna_ch.multiMap { meta, seqs, db ->
        seqs: [meta, seqs]
        db: db
    }
    // rrna_ch.view{ "rrna_ch - ${it}" }
    // ch_test = rrna_chs.seqs.map { meta, reads ->
    //     [meta + ['rrna_read_count': reads.countFasta()], reads]
    // }
    // ch_test.view{ meta, _reads -> "ch_test - [${meta.id}, ${meta.platform}, ${meta.single_end}] - ${meta.rrna_read_count}" }

    MAPSEQ_OTU_KRONA(rrna_chs.seqs, rrna_chs.db)
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)


    // Pfam profiling
    pfam_db = dbs.pfam
        .map { meta, fp ->
            file("${fp}/${meta.base_dir}/${meta.files.hmm}")
        }
        .first()

    PROFILE_HMMSEARCH_PFAM(
        READSMERGE.out.reads_fasta,
        pfam_db
    )
    ch_versions = ch_versions.mix(PROFILE_HMMSEARCH_PFAM.out.versions)


    // MultiQC
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    trim_meta = { meta, v -> [[meta.id, meta.single_end, meta.platform], v] }
    // decontam_stats.map(trim_meta).view{ "decontam_stats - ${it}" }
    // QC.out.fastp_json.map(trim_meta).view{ "QC.out.fastp_json - ${it}" }
    multiqc_ch = qc_stats.map(trim_meta)
        .join(decontam_stats.map(trim_meta), remainder: true)
        .multiMap { it ->
            names: it[0][0]
            files: it[1..-1]
        }
    // multiqc_ch.names.view { "multiqc_ch.names - ${it}" }

    MULTIQC(
        multiqc_ch.files.flatten().filter{ it }.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    // Collate software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'ASSEMBLY_ANALYSIS_PIPELINE_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { collated_versions }

    } // end download_dbs condition

    emit:
    versions = ch_versions  // channel: [ path(versions.yml) ]
}
