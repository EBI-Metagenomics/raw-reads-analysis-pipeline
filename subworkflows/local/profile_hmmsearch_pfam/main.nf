include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/seqkit/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { SEEDED_HMMER } from '../seeded_hmmer/main'
include { FRAGGENESCANRS } from '../../../modules/local/fraggenescanrs/main'
include { SELECTFRAMES } from '../../../modules/local/selectframes/main'

workflow PROFILE_HMMSEARCH_PFAM {
    take:
    reads_fasta
    pfam_db

    main:
    ch_versions = Channel.empty()

    // The seeded engines search a whole HMM file, so they have no Pfam chunking to spread a
    // sample over: without this, translation, seeding and HMMER are each one task per sample over
    // the whole library. Chunk the reads instead, before translation, so ORF counts per task stay
    // bounded whatever the library's depth. hmmsearch keeps chunking its ORFs by size below, since
    // it is already spread over the database's chunks.
    ch_reads = reads_fasta
    if (params.pfam_search_engine in ['nail', 'kfpfam'] && params.pfam_seeded_chunk_reads > 0) {
        ch_reads = reads_fasta
            .splitFasta(by: params.pfam_seeded_chunk_reads, elem: 1, file: true)
            .groupTuple()
            .flatMap { meta, fastas ->
                def chunks = fastas instanceof Collection ? fastas : [fastas]
                // read_chunk makes each task's key unique; read_chunks is how many to group back
                chunks.withIndex().collect { fasta, i ->
                    tuple(meta + [read_chunk: i + 1, read_chunks: chunks.size()], fasta)
                }
            }
    }

    SEQKIT_TRANSLATE(ch_reads)

    ch_orfs = SEQKIT_TRANSLATE.out.fastx

    if (params.pfam_orf_filter) {
        // Reads with no sequence have no ORFs either, and FragGeneScanRs needs input
        FRAGGENESCANRS(ch_reads.filter { _meta, fasta -> fasta.countFasta(limit: 1) > 0 })
        SELECTFRAMES(ch_reads.join(SEQKIT_TRANSLATE.out.fastx).join(FRAGGENESCANRS.out.calls))
        ch_versions = ch_versions.mix(FRAGGENESCANRS.out.versions, SELECTFRAMES.out.versions)
        ch_orfs = SELECTFRAMES.out.orfs
    }

    ch_alignments = Channel.empty()

    if (params.pfam_search_engine in ['nail', 'kfpfam']) {
        if (params.pfam_search_engine == 'kfpfam' && !(params.kfpfam_prep && params.kfpfam_index)) {
            error(
                "pfam_search_engine 'kfpfam' needs kfpfam_prep and kfpfam_index, the prep directory " +
                "and the tiled k-mer index built from the same Pfam-A.hmm this run searches. " +
                "hpc/pfam-arms/kfpfam_setup.sh builds both, and kfpfam_sif the image that runs them."
            )
        }

        SEEDED_HMMER(ch_orfs, pfam_db)
        ch_versions = ch_versions.mix(SEEDED_HMMER.out.versions)
        ch_alignments = SEEDED_HMMER.out.alignments
        ch_domtbl = SEEDED_HMMER.out.domtbl
    }
    else {
        // chunk seqs
        ch_chunked_fasta = ch_orfs
            .splitFasta(
                size: params.hmmsearch_chunksize,
                elem: 1,
                file: true
            )

        // combine after chunking to allow for a chunked pfam database
        ch_chunked_pfam_in = ch_chunked_fasta
            .combine(pfam_db)
            .map { meta, reads, db -> tuple(meta, tuple(reads, db)) }
    
        ch_chunked_pfam_in = ch_chunked_pfam_in
            .groupTuple()  // group to count the number of chunks
            .flatMap {  // un-group with information about the number of chunks (to aid re-grouping)
                meta, chunks ->
                def chunks_ = chunks instanceof Collection ? chunks : [chunks]
                def chunksize = chunks_.size()
                return chunks_.collect {
                    chunk ->
                    tuple(groupKey(meta, chunksize), chunk)
                }
            }
            .map { meta, v ->
                def (reads, db) = v
                return [meta, db, reads, false, true, true] 
            }

        HMMER_HMMSEARCH(ch_chunked_pfam_in)
        ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

        COMBINEHMMSEARCHTBL(
            HMMER_HMMSEARCH.out.domain_summary.groupTuple()
        )
        ch_domtbl = COMBINEHMMSEARCHTBL.out.concatenated_result
    }

    PARSEHMMSEARCHCOVERAGE(ch_domtbl)
    ch_versions = ch_versions.mix(PARSEHMMSEARCHCOVERAGE.out.versions)

    emit:
    profile = PARSEHMMSEARCHCOVERAGE.out.tsv
    stats = PARSEHMMSEARCHCOVERAGE.out.stats
    alignments = ch_alignments
    versions = ch_versions // channel: [ versions.yml ]
}
