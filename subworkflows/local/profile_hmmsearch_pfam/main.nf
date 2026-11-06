include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/seqkit/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { NAIL_HMMER } from '../nail_hmmer/main'

workflow PROFILE_HMMSEARCH_PFAM {
    take:
    reads_fasta
    pfam_db

    main:
    ch_versions = Channel.empty()

    SEQKIT_TRANSLATE(reads_fasta)

    ch_alignments = Channel.empty()

    if (params.pfam_search_engine == 'nail') {
        NAIL_HMMER(SEQKIT_TRANSLATE.out.fastx, pfam_db)
        ch_versions = ch_versions.mix(NAIL_HMMER.out.versions)
        ch_alignments = NAIL_HMMER.out.alignments
        ch_domtbl = NAIL_HMMER.out.domtbl
    }
    else {
        // chunk seqs
        ch_chunked_fasta = SEQKIT_TRANSLATE.out.fastx
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
