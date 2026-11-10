include { NAIL_SEED } from '../../../modules/local/nail_seed/main'
include { HMMER_PAIRS } from '../../../modules/local/hmmer_pairs/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { MERGEPARQUET } from '../../../modules/local/mergeparquet/main'

workflow NAIL_HMMER {
    take:
    orfs // channel: [ val(meta), path(fasta) ]
    pfam_db // channel: path(hmm), one per chunk of a chunked Pfam, else a single element

    main:
    // nail stops on an empty target file; hmmsearch's chunking also skips samples with no ORFs
    ch_orfs = orfs.filter { _meta, fasta -> fasta.countFasta(limit: 1) > 0 }

    // nail searches one whole HMM file, so a chunked Pfam means one search per chunk. Sorted
    // by path so the concatenated domtbl doesn't depend on the order the chunks arrive in.
    ch_hmm_chunks = pfam_db
        .toSortedList()
        .flatMap { hmms -> hmms.withIndex().collect { hmm, i -> tuple(i + 1, hmms.size(), hmm) } }

    // hmm_chunk makes the key unique per (sample, chunk), so the seeds join back to their own
    // chunk; hmm_chunks is how many pieces to group back into one profile.
    ch_chunks = ch_orfs
        .combine(ch_hmm_chunks)
        .map { meta, fasta, chunk, chunks, hmm ->
            tuple(meta + [hmm_chunk: chunk, hmm_chunks: chunks], fasta, hmm)
        }

    NAIL_SEED(ch_chunks)

    // num_models is the whole database either way, so E-values match an unchunked search
    HMMER_PAIRS(
        ch_chunks.join(NAIL_SEED.out.seeds),
        params.databases.pfam.variables.num_models,
    )

    // Back to one domtbl and one Parquet per sample. Both steps also run for a single chunk,
    // so the chunked and unchunked paths stay the same shape.
    ch_grouped = HMMER_PAIRS.out.domtbl
        .join(HMMER_PAIRS.out.alignments)
        .map { meta, domtbl, parquet ->
            def key = meta.findAll { k, _v -> !(k in ['hmm_chunk', 'hmm_chunks']) }
            tuple(groupKey(key, meta.hmm_chunks), domtbl, parquet)
        }
        .groupTuple()

    COMBINEHMMSEARCHTBL(ch_grouped.map { meta, domtbls, _parquets -> tuple(meta, domtbls) })
    MERGEPARQUET(ch_grouped.map { meta, _domtbls, parquets -> tuple(meta, parquets) })

    emit:
    domtbl = COMBINEHMMSEARCHTBL.out.concatenated_result // channel: [ val(meta), path(domtbl.gz) ]
    alignments = MERGEPARQUET.out.parquet // channel: [ val(meta), path(parquet) ]
    versions = NAIL_SEED.out.versions.mix(HMMER_PAIRS.out.versions, MERGEPARQUET.out.versions) // channel: [ versions.yml ]
}
