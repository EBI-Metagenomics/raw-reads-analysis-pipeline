include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/seqkit/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { NAIL_HMMER } from '../nail_hmmer/main'
include { FRAGGENESCANRS } from '../../../modules/local/fraggenescanrs/main'
include { SELECTFRAMES } from '../../../modules/local/selectframes/main'

workflow PROFILE_HMMSEARCH_PFAM {
    take:
    reads_fasta
    pfam_db

    main:
    ch_versions = Channel.empty()

    SEQKIT_TRANSLATE(reads_fasta)

    ch_orfs = SEQKIT_TRANSLATE.out.fastx

    if (params.pfam_orf_filter) {
        // Reads with no sequence have no ORFs either, and FragGeneScanRs needs input
        FRAGGENESCANRS(reads_fasta.filter { _meta, fasta -> fasta.countFasta(limit: 1) > 0 })
        SELECTFRAMES(reads_fasta.join(SEQKIT_TRANSLATE.out.fastx).join(FRAGGENESCANRS.out.calls))
        ch_versions = ch_versions.mix(FRAGGENESCANRS.out.versions, SELECTFRAMES.out.versions)
        ch_orfs = SELECTFRAMES.out.orfs
    }

    ch_alignments = Channel.empty()

    if (params.pfam_search_engine == 'nail') {
        // nail searches one whole HMM file: a chunked Pfam would give it one chunk and a
        // partial profile with no error. Fail here rather than emitting a wrong answer.
        if (params.databases.pfam.containsKey('chunked') && params.databases.pfam.chunked) {
            error(
                "pfam_search_engine 'nail' needs an unchunked Pfam database, but " +
                "params.databases.pfam.chunked is set. Point params.databases.pfam at a single " +
                "Pfam-A.hmm (remote_path or local_path plus files.hmm, no chunked = true), or " +
                "use pfam_search_engine 'hmmsearch', which searches each chunk."
            )
        }

        NAIL_HMMER(ch_orfs, pfam_db)
        ch_versions = ch_versions.mix(NAIL_HMMER.out.versions)
        ch_alignments = NAIL_HMMER.out.alignments
        ch_domtbl = NAIL_HMMER.out.domtbl
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
