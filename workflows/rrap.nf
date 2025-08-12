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
include { ADDHEADER as ADDHEADER_RRNA } from '../modules/local/addheader/main'
include { ADDHEADER as ADDHEADER_MOTUS } from '../modules/local/addheader/main'
include { BBMAP_REFORMAT_STANDARDISE } from '../modules/local/bbmap/reformat_standardise/main'
include { SEQKIT_SHUFFLE_FASTA } from '../modules/local/seqkit_shuffle_fasta/main'

include { RRNA_EXTRACTION } from '../subworkflows/local/rrna_extraction/main'
include { MAPSEQ_OTU_KRONA } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MULTIQC as MULTIQC_RUN } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_STUDY } from '../modules/nf-core/multiqc/main'
include { PROFILE_HMMSEARCH_PFAM } from '../subworkflows/local/profile_hmmsearch_pfam/main'
include { samplesheetToList } from 'plugin/nf-schema'


workflow PIPELINE {
    main:
    println(new File("${projectDir}/assets/header.txt").text)

    ch_versions = Channel.empty()

   // Fetch databases
    db_ch = Channel
        .from(
            params.databases.collect { k, v ->
                if (v instanceof Map) {
                    if (v.containsKey('chunked') && v['chunked']) {
                        v.collect { k_, v_ ->
                            if (v_ instanceof Map) {
                                if (v_.containsKey('files')) {
                                    return [id: k, chunk_id: k_] + v_
                                }
                            }
                        }
                    } else if (v.containsKey('files')) {
                        return [id: k] + v
                    }
                }
            }.flatten()
        )
        .filter { it }

    FETCHDB(db_ch, "${launchDir}/${params.databases.cache_path}")
    dbs_path_ch = FETCHDB.out.dbs

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

    // Parse samplesheet and fetch reads
    def groupReads = { sample, fq1, fq2, single_end, instrument_platform ->
        def single_file = (fq2 == [])
        [
            [
                'id': sample,
                'single_end': single_end,
                'interleaved': (!single_end) && single_file,
                'instrument_platform': instrument_platform,
            ],
            single_file ? [file(fq1)] : [file(fq1), file(fq2)],
        ]
    }
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))

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

    classified_reads
    .filter { meta, _reads -> (meta.long_reads && (!meta.single_end)) }
    .collect { meta, _reads ->
        error "Error: Long reads (ONT or PB) cannot be paired-end — ${meta}"
    }

    // De-interleave interleaved paired-end reads
    BBMAP_REFORMAT_STANDARDISE(classified_reads, 'fastq.gz')
    classified_reads = BBMAP_REFORMAT_STANDARDISE.out.reformated

    // QC
    if (params.skip_qc) {
        classified_reads.set { qc_reads }
        qc_stats = Channel.empty()
    }
    else {
        QC(classified_reads)
        ch_versions = ch_versions.mix(QC.out.versions)

        qc_reads = QC.out.fastq
        qc_stats = QC.out.fastp_json

        qc_read_counts = qc_stats.map {
            meta, json_file ->
            def json = new groovy.json.JsonSlurper().parseText(json_file.text)
            return tuple(
                meta,
                tuple(
                    json["summary"]["before_filtering"]["total_reads"],
                    json["summary"]["after_filtering"]["total_reads"],
                )
            )
        }

        qc_reads = qc_reads.join(qc_read_counts)
        .map{ meta, reads, counts ->
            tuple(
                meta + ['read_count': counts[0], 'qc_read_count': counts[1]],
                reads
            )
        }
        .filter{ meta, _reads -> meta.qc_read_count > 0 }
    }


    // DECONTAMINATION
    decontam_stats = Channel.empty()

    if (params.skip_decontam) {
        qc_reads.set { clean_reads }
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

        clean_reads = DECONTAM_SHORTREAD.out.decontaminated_reads.mix(DECONTAM_LONGREAD.out.decontaminated_reads)

    }

    // Merge reads
    READSMERGE(clean_reads)
    merged_reads = READSMERGE.out.reads_fasta

    // Get post-decontam stats
    decontam_stats = READSMERGE.out.fastp_summary_json

    decontam_read_counts = decontam_stats.map {
        meta, json_file ->
        def json = new groovy.json.JsonSlurper().parseText(json_file.text)
        return tuple(
            meta,
            tuple(
                json["summary"]["before_filtering"]["total_reads"],
                json["summary"]["after_filtering"]["total_reads"],
            )
        )
    }

    clean_reads = clean_reads
        .join(decontam_read_counts)
        .map{ meta, reads, counts ->
            tuple(
                meta + [
                    'clean_read_count': counts[0],
                    'merged_read_count': counts[1],
                ],
                reads
            )
        }
    .filter{ meta, _reads -> meta.clean_read_count > 0 }

    merged_reads = merged_reads
        .join(decontam_read_counts)
        .map{ meta, reads, counts ->
            tuple(
                meta + [
                    'clean_read_count': counts[0],
                    'merged_read_count': counts[1],
                ],
                reads
            )
        }
        .filter{ meta, _reads -> meta.merged_read_count > 0 }

    // mOTUs
    motus_db = dbs.motus
        .map { meta, fp ->
            file("${fp}/${meta.files.base_dir}")
        }
        .first()

    MOTUS_KRONA(
        clean_reads,
        motus_db,
    )
    ch_versions = ch_versions.mix(MOTUS_KRONA.out.versions)

    ADDHEADER_MOTUS(
        MOTUS_KRONA.out.krona,
        "# ${params.results_file_headers.motus_taxonomy.join('\t')}",
    )


    // rrna_extraction
    rfam_db = dbs.rfam
        .map { meta, fp ->
            file("${fp}/${meta.files.models}")
        }
        .first()

    claninfo_db = dbs.rfam
        .map { meta, fp ->
            file("${fp}/${meta.files.claninfo}")
        }
        .first()

    RRNA_EXTRACTION(
        merged_reads,
        rfam_db,
        claninfo_db,
    )
    ch_versions = ch_versions.mix(RRNA_EXTRACTION.out.versions)

    lsu_db = dbs.silva_lsu
        .map { meta, fp ->
            [
                [
                    file("${fp}/${meta.files.fasta}"),
                    file("${fp}/${meta.files.tax}"),
                    file("${fp}/${meta.files.otu}"),
                    file("${fp}/${meta.files.mscluster}"),
                    meta.id,
                ]
            ]
        }
        .first()

    ssu_db = dbs.silva_ssu
        .map { meta, fp ->
            [
                [
                    file("${fp}/${meta.files.fasta}"),
                    file("${fp}/${meta.files.tax}"),
                    file("${fp}/${meta.files.otu}"),
                    file("${fp}/${meta.files.mscluster}"),
                    meta.id,
                ]
            ]
        }
        .first()

    lsu_ch = RRNA_EXTRACTION.out.lsu_fasta
        .map { meta, fp -> [meta + ['db_label': 'SILVA-LSU'], fp] }
        .combine(lsu_db)
    ssu_ch = RRNA_EXTRACTION.out.ssu_fasta
        .map { meta, fp -> [meta + ['db_label': 'SILVA-SSU'], fp] }
        .combine(ssu_db)
    rrna_ch = lsu_ch.mix(ssu_ch)
    rrna_chs = rrna_ch.multiMap { meta, seqs, db ->
        seqs: [meta, seqs]
        db: db
    }

    MAPSEQ_OTU_KRONA(rrna_chs.seqs, rrna_chs.db)
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)

    ADDHEADER_RRNA(
        MAPSEQ_OTU_KRONA.out.krona_input,
        "# ${params.results_file_headers.silva_taxonomy.join('\t')}",
    )


    if (params.skip_functional) {
        pfam_status = merged_reads.map { meta, _fp -> [meta.id, false] }
    }
    else {
        // Pfam profiling
        pfam_db = dbs.pfam
            .map { meta, fp ->
                file("${fp}/${meta.files.hmm}")
            }

        if (params.hmmsearch_subsampling == -1) {
            pfam_reads = merged_reads
        }
        else {
            SEQKIT_SHUFFLE_FASTA(merged_reads, params.hmmsearch_subsampling, true)
            pfam_reads = SEQKIT_SHUFFLE_FASTA.out.fasta
        }

        PROFILE_HMMSEARCH_PFAM(
            pfam_reads,
            pfam_db,
        )
        ch_versions = ch_versions.mix(PROFILE_HMMSEARCH_PFAM.out.versions)
        
        pfam_status = PROFILE_HMMSEARCH_PFAM.out.profile.map { meta, fp -> [meta.id, fp.exists() && (fp.readLines().size() > 0)] }
    }

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

    trim_meta = { meta, v ->
        [
            [
                id: meta.id,
                single_end: meta.single_end,
                instrument_platform: meta.instrument_platform,
            ],
            v,
        ]
    }

    // per Run
    multiqc_run_ch = qc_stats
        .map(trim_meta)
        .mix(decontam_stats.map(trim_meta))
        .groupTuple()

    MULTIQC_RUN(
        multiqc_run_ch,
        [],
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    // Study
    multiqc_study_ch = qc_stats
        .map(trim_meta)
        .mix(decontam_stats.map(trim_meta))
        .groupTuple()
        .multiMap { meta, files ->
            names: (1..files.size()).collect { meta.id }
            files: files
        }
    multiqc_study_ch = Channel
        .value([id: "study"])
        .combine(multiqc_study_ch.files.flatten().collect())
        .map { new Tuple(it[0], it[1..-1]) }

    MULTIQC_STUDY(
        multiqc_study_ch,
        [],
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
            name: 'RAW_READS_ANALYSIS_PIPELINE_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { collated_versions }

    reads_status = classified_reads.map { meta, _reads -> [meta.id, meta.read_count > 0] }

    qc_status = qc_reads.map { meta, _reads -> [meta.id, meta.qc_read_count > 0] }

    decontam_status = clean_reads.map { meta, _reads -> [meta.id, meta.clean_read_count > 0] }

    motus_status = ADDHEADER_MOTUS.out.file_with_header.map { meta, fp -> [meta.id, fp.exists() && (fp.readLines().size() > 0)] }

    silvassu_status = ADDHEADER_RRNA.out.file_with_header
        .filter { meta, _fp -> meta.db_label == 'SILVA-SSU' }
        .map { meta, fp -> [meta.id, fp.exists() && (fp.readLines().size() > 0)] }

    silvalsu_status = ADDHEADER_RRNA.out.file_with_header
        .filter { meta, _fp -> meta.db_label == 'SILVA-LSU' }
        .map { meta, fp -> [meta.id, fp.exists() && (fp.readLines().size() > 0)] }

    run_status = reads_status
        .join(qc_status, remainder: true)
        .join(decontam_status, remainder: true)
        .join(motus_status, remainder: true)
        .join(silvassu_status, remainder: true)
        .join(silvalsu_status, remainder: true)
        .join(pfam_status, remainder: true)

    run_status
        .filter { _meta_id, _reads, qc, _decontam, _motus, _silvassu, _silvalsu, _pfam -> qc }
        .map { meta_id, _reads, _qc, decontam, motus, silvassu, silvalsu, pfam ->
            {
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
            }
        }
        .collectFile(name: "qc_passed_runs.csv", storeDir: params.outdir, newLine: true, cache: false)

    run_status
        .filter { _meta_id, _reads, qc, _decontam, _motus, _silvassu, _silvalsu, _pfam -> !qc }
        .map { meta_id, _reads, _qc, decontam, motus, silvassu, silvalsu, pfam ->
            {
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
            }
        }
        .collectFile(name: "qc_failed_runs.csv", storeDir: params.outdir, newLine: true, cache: false)

    run_status
        .map { meta_id, reads, qc, decontam, motus, silvassu, silvalsu, pfam ->
            {
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
                return "${meta_id},${status},${reads ? "reads_yes" : "read_no"},${qc ? "qc_yes" : "qc_no"},${decontam ? "decontam_yes" : "decontam_no"},${motus ? "motus_yes" : "motus_no"},${silvassu ? "silva-ssu_yes" : "silva-ssu_no"},${silvalsu ? "silva-lsu_yes" : "silva-lsu_no"},${pfam ? "pfam_yes" : "pfam_no"}"
            }
        }
        .collectFile(name: "run_status.csv", storeDir: params.outdir, newLine: true, cache: false)

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
    collated_versions = collated_versions // channel: [ path(versions.yml) ]
}
