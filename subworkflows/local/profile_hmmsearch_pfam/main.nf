include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/seqkit/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/local/hmmer/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'

workflow PROFILE_HMMSEARCH_PFAM {
    take:
    reads_fasta
    pfam_db

    main:
    ch_versions = Channel.empty()

    SEQKIT_TRANSLATE(reads_fasta)

    ch_chunked_fasta = SEQKIT_TRANSLATE.out.fastx
        .splitFasta(
            size: params.hmmsearch_chunksize,
            elem: 1,
            file: true
        )
        .groupTuple()
        .flatMap {
            meta, chunks ->
            def chunks_ = chunks instanceof Collection ? chunks : [chunks]
            def chunksize = chunks_.size()
            return chunks_.collect {
                chunk ->
                tuple(groupKey(meta, chunksize), chunk)
            }
        }

    ch_chunked_pfam_in = ch_chunked_fasta
        .combine(pfam_db)
        .map { meta, reads, db -> [meta, db, reads, false, true, true] }

    HMMER_HMMSEARCH(ch_chunked_pfam_in)
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    COMBINEHMMSEARCHTBL(
        HMMER_HMMSEARCH.out.domain_summary.groupTuple()
    )

    COMBINEHMMSEARCHTBL.out.concatenated_result.view { "hmmsearch_concat - ${it}" }

    PARSEHMMSEARCHCOVERAGE(COMBINEHMMSEARCHTBL.out.concatenated_result)
    ch_versions = ch_versions.mix(PARSEHMMSEARCHCOVERAGE.out.versions)

    emit:
    profile = PARSEHMMSEARCHCOVERAGE.out.tsv
    stats = PARSEHMMSEARCHCOVERAGE.out.stats
    versions = ch_versions // channel: [ versions.yml ]
}
