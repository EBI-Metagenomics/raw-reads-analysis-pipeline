include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/seqkit/translate/main'
include { FASTAEMBEDLENGTH } from '../../../modules/local/fastaembedlength/main'
include { HMMER_HMMSEARCH } from '../../../modules/local/hmmer/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { CHUNKFASTX} from  '../../../modules/local/chunkfastx/main'

workflow PROFILE_HMMSEARCH_PFAM {
    take:
    reads_fasta
    pfam_db

    main:
    ch_versions = Channel.empty()

    FASTAEMBEDLENGTH(reads_fasta)

    SEQKIT_TRANSLATE(FASTAEMBEDLENGTH.out.fasta)

    CHUNKFASTX(SEQKIT_TRANSLATE.out.fastx)

    ch_chunked_fasta = CHUNKFASTX.out.chunked_reads.flatMap {
        meta, chunks ->
        def chunks_ = chunks instanceof Collection ? chunks : [chunks]
        return chunks_.collect {
            chunk ->
            tuple(groupKey(meta, chunks_.size()), chunk)
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

    PARSEHMMSEARCHCOVERAGE(COMBINEHMMSEARCHTBL.out.concatenated_result)
    ch_versions = ch_versions.mix(PARSEHMMSEARCHCOVERAGE.out.versions)

    emit:
    profile = PARSEHMMSEARCHCOVERAGE.out.tsv
    versions = ch_versions // channel: [ versions.yml ]
}
