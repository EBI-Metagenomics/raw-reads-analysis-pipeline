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
include { ADDHEADER as ADDHEADER_RRNA } from '../modules/local/addheader/main'
include { ADDHEADER as ADDHEADER_MOTUS } from '../modules/local/addheader/main'

include { RRNA_EXTRACTION } from '../subworkflows/local/rrna_extraction/main'
include { MAPSEQ_OTU_KRONA } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { MULTIQC as MULTIQC_RUN } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_STUDY } from '../modules/nf-core/multiqc/main'
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
    FETCHDB(db_ch, "${projectDir}/${params.databases.cache_path}")
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
    def groupReads = { study, sample, fq1, fq2, library_layout, library_strategy, instrument_platform ->
        [
            [
                'id': sample,
                'study': study,
                'single_end': fq2 == [] ? true : false,
                'library_layout': library_layout,
                'library_strategy': library_strategy,
                'instrument_platform': instrument_platform
            ],
            fq2 == [] ? file(fq1) : [file(fq1), file(fq2)]
        ]
    }
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))

    // [ study, sample, read1, [read2], library_layout, library_strategy, instrument_platform]
    fetch_reads_transformed = samplesheet.map(groupReads)

    classified_reads = fetch_reads_transformed.map { meta, reads ->
        // Long reads
        if (["ONT", "PB"].contains(meta.instrument_platform)) {
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

    // Get read count per fastq row
    clean_reads = clean_reads.map { meta, reads ->
        [meta + ['clean_read_count': (meta.single_end ? reads : reads[0]).countFastq()], reads]
    }
    // clean_reads.view{ meta, _reads -> "clean_reads - [${meta.id}, ${meta.instrument_platform}, ${meta.single_end}] - ${meta.clean_read_count}" }
    clean_reads = clean_reads.filter { meta, _reads ->
        meta.clean_read_count > 0
    }

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

    ADDHEADER_MOTUS(
        MOTUS_KRONA.out.krona,
        "# ${params.results_file_headers.motus_taxonomy.join('\t')}"
    )

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

    lsu_ch = RRNA_EXTRACTION.out.lsu_fasta
        .map { meta, fp -> [meta+['db_label': 'SILVA-LSU'], fp]}
        .combine(lsu_db)
    ssu_ch = RRNA_EXTRACTION.out.ssu_fasta
        .map { meta, fp -> [meta+['db_label': 'SILVA-SSU'], fp]}
        .combine(ssu_db)
    rrna_ch = lsu_ch.mix(ssu_ch)
    rrna_chs = rrna_ch.multiMap { meta, seqs, db ->
        seqs: [meta, seqs]
        db: db
    }
    // rrna_ch.view{ "rrna_ch - ${it}" }
    // ch_test = rrna_chs.seqs.map { meta, reads ->
    //     [meta + ['rrna_read_count': reads.countFasta()], reads]
    // }
    // ch_test.view{ meta, _reads -> "ch_test - [${meta.id}, ${meta.instrument_platform}, ${meta.single_end}] - ${meta.rrna_read_count}" }

    MAPSEQ_OTU_KRONA(rrna_chs.seqs, rrna_chs.db)
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)

    ADDHEADER_RRNA(
        MAPSEQ_OTU_KRONA.out.krona_input,
        "# ${params.results_file_headers.silva_taxonomy.join('\t')}"
    )


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

    trim_meta = { meta, v -> [[meta.id, meta.single_end, meta.instrument_platform], v] }
    // decontam_stats.map(trim_meta).view{ "decontam_stats - ${it}" }
    // QC.out.fastp_json.map(trim_meta).view{ "QC.out.fastp_json - ${it}" }
    multiqc_ch = qc_stats.map(trim_meta)
        .join(decontam_stats.map(trim_meta), remainder: true)
    // multiqc_ch.view { "multiqc_ch - ${it}" }

    // per Run
    // MULTIQC_RUN(
    //     multiqc_ch,
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     [],
    // )

    // Study
    multiqc_study_ch = multiqc_ch
        .multiMap { it ->
            names: it[0][0]
            files: it[1..-1]
        }
    // multiqc_study_ch.names.view { "multiqc_study_ch.names - ${it}" }

    MULTIQC_STUDY(
        multiqc_study_ch.files.flatten().filter{ it }.collect(),
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

    reads_status = classified_reads
        .map{ meta, _reads -> [meta.id, meta.read_count > 0] }
    // reads_status.view { "reads_status - ${it}" }

    qc_status = qc_reads
        .map{ meta, _reads -> [meta.id, meta.qc_read_count > 0] }
    // qc_status.view { "qc_status - ${it}" }

    decontam_status = clean_reads
        .map{ meta, _reads -> [meta.id, meta.clean_read_count > 0] }
    // decontam_status.view { "decontam_status - ${it}" }

    motus_status = ADDHEADER_MOTUS.out.file_with_header
        .map{ meta, fp -> [meta.id, fp.exists() && (fp.readLines().size()>0)] }
    // motus_status.view { "motus_status - ${it}" }

    silvassu_status = ADDHEADER_RRNA.out.file_with_header
        .filter{ meta, _fp -> meta.db_label=='SILVA-SSU' }
	    .map{ meta, fp -> [meta.id, fp.exists() && (fp.readLines().size()>0)] }
    // silvassu_status.view { "silvassu_status - ${it}" }

    silvalsu_status = ADDHEADER_RRNA.out.file_with_header
        .filter{ meta, _fp -> meta.db_label=='SILVA-LSU' }
        .map{ meta, fp -> [meta.id, fp.exists() && (fp.readLines().size()>0)] }
    // silvalsu_status.view { "silvalsu_status - ${it}" }

    pfam_status = PROFILE_HMMSEARCH_PFAM.out.profile
        .map{ meta, fp -> [meta.id, fp.exists() && (fp.readLines().size()>0)] }
    // pfam_status.view { "pfam_status - ${it}" }

    run_status = reads_status
        .join( qc_status, remainder: true )
        .join( decontam_status, remainder: true )
        .join( motus_status, remainder: true )
        .join( silvassu_status, remainder: true )
        .join( silvalsu_status, remainder: true )
        .join( pfam_status, remainder: true )
        .map{ meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam -> {
	        def status = "all_results"
	        if (decontam == false) {
	            status = "no_reads"
	        }
	        if (![motus, silvassu, silvalsu, pfam].any()) {
	            status = "no_results"
	        }
	        if (![motus, silvassu, silvalsu, pfam].every()) {
	            status = "missing_results"
	        }
	        return "${meta_id},${status},${reads ? "reads_yes":"read_no"},${qc ? "qc_yes":"qc_no"},${decontam ? "decontam_yes":"decontam_no"},${motus ? "motus_yes":"motus_no"},${silvassu ? "silva-ssu_yes":"silva-ssu_no"},${silvalsu ? "silva-lsu_yes":"silva-lsu_no"},${pfam ? "pfam_yes":"pfam_no"}"
	    }}

    run_status
        .filter{ meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam -> qc }
        .map{ meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam -> {
	        def status = "all_results"
	        if (decontam == false) {
	            status = "no_reads"
	        }
	        if (![motus, silvassu, silvalsu, pfam].any()) {
	            status = "no_results"
	        }
	        if (![motus, silvassu, silvalsu, pfam].every()) {
	            status = "missing_results"
	        }
	        return "${meta_id},${status}"
	    }}
        .collectFile(name: "qc_passed_runs.csv", storeDir: params.outdir, newLine: true, cache: false)

    run_status
        .filter{ meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam -> !qc }
        .map{ meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam -> {
	        def status = "all_results"
	        if (decontam == false) {
	            status = "no_reads"
	        }
	        if (![motus, silvassu, silvalsu, pfam].any()) {
	            status = "no_results"
	        }
	        if (![motus, silvassu, silvalsu, pfam].every()) {
	            status = "missing_results"
	        }
	        return "${meta_id},${status}"
	    }}
        .collectFile(name: "qc_failed_runs.csv", storeDir: params.outdir, newLine: true, cache: false)

    run_status
        .map{ meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam -> {
	        def status = "all_results"
	        if (decontam == false) {
	            status = "no_reads"
	        }
	        if (![motus, silvassu, silvalsu, pfam].any()) {
	            status = "no_results"
	        }
	        if (![motus, silvassu, silvalsu, pfam].every()) {
	            status = "missing_results"
	        }
	        return "${meta_id},${status},${reads ? "reads_yes":"read_no"},${qc ? "qc_yes":"qc_no"},${decontam ? "decontam_yes":"decontam_no"},${motus ? "motus_yes":"motus_no"},${silvassu ? "silva-ssu_yes":"silva-ssu_no"},${silvalsu ? "silva-lsu_yes":"silva-lsu_no"},${pfam ? "pfam_yes":"pfam_no"}"
	    }}
        .collectFile(name: "run_status.csv", storeDir: params.outdir, newLine: true, cache: false)

    emit:
    versions = ch_versions                                // channel: [ path(versions.yml) ]
    pfam_profile = PROFILE_HMMSEARCH_PFAM.out.profile     // channel: [ meta, path ]
    rrna_profile = ADDHEADER_RRNA.out.file_with_header    // channel: [ meta, path ]
    motus_profile = ADDHEADER_MOTUS.out.file_with_header  // channel: [ meta, path ]
    decontam_stats = decontam_stats                       // channel: [ meta, path ]
    qc_stats = qc_stats                                   // channel: [ meta, path ]
}
