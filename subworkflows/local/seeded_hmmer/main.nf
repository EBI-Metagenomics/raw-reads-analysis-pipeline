include { NAIL_SEED } from '../../../modules/local/nail_seed/main'
include { KFPFAM_FILTER } from '../../../modules/local/kfpfam_filter/main'
include { HMMER_PAIRS } from '../../../modules/local/hmmer_pairs/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { MERGEPARQUET } from '../../../modules/local/mergeparquet/main'

workflow SEEDED_HMMER {
    take:
    orfs // channel: [ val(meta), path(fasta) ]
    pfam_db // channel: path(hmm), one per chunk of a chunked Pfam, else a single element

    main:
    // Both engines stop on an empty target file; hmmsearch's chunking also skips samples with no ORFs
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

    ch_versions = Channel.empty()

    if (params.pfam_search_engine == 'kfpfam') {
        // kfpfam's index covers the whole database, so it seeds once per sample whatever the Pfam
        // chunking is. Its seeds name models from every chunk; a chunk's HMMER_PAIRS task reads the
        // HMM file it was given and ignores the rest.
        KFPFAM_FILTER(
            ch_orfs,
            file(params.kfpfam_prep, checkIfExists: true),
            file(params.kfpfam_index, checkIfExists: true),
            params.kfpfam_bar,
            params.kfpfam_gate,
            params.kfpfam_min_gapped,
            params.kfpfam_block,
            params.databases.pfam.variables.num_models,
        )
        ch_versions = ch_versions.mix(KFPFAM_FILTER.out.versions)
        // One seeds file per read chunk, copied onto that chunk's Pfam-chunk keys. The key is
        // (sample, read chunk), since PROFILE_HMMSEARCH_PFAM may have split the reads.
        ch_seeds = ch_chunks
            .map { meta, _fasta, _hmm -> tuple([meta.id, meta.read_chunk], meta) }
            .combine(
                KFPFAM_FILTER.out.seeds.map { meta, seeds -> tuple([meta.id, meta.read_chunk], seeds) },
                by: 0,
            )
            .map { _key, meta, seeds -> tuple(meta, seeds) }
    }
    else {
        // nail's seeding reads one HMM file, so it runs per (sample, chunk)
        NAIL_SEED(ch_chunks)
        ch_versions = ch_versions.mix(NAIL_SEED.out.versions)
        ch_seeds = NAIL_SEED.out.seeds
    }

    // num_models is the whole database either way, so E-values match an unchunked search
    HMMER_PAIRS(
        ch_chunks.join(ch_seeds),
        params.databases.pfam.variables.num_models,
    )
    ch_versions = ch_versions.mix(HMMER_PAIRS.out.versions)

    // Back to one domtbl and one Parquet per sample: a task per (read chunk, Pfam chunk), so the
    // group is as big as their product. Both steps also run for a single chunk of each, so the
    // chunked and unchunked paths stay the same shape. Sorted so the concatenated domtbl does not
    // depend on the order the chunks finish in.
    ch_chunk_keys = ['hmm_chunk', 'hmm_chunks', 'read_chunk', 'read_chunks']
    ch_grouped = HMMER_PAIRS.out.domtbl
        .join(HMMER_PAIRS.out.alignments)
        .map { meta, domtbl, parquet ->
            def key = meta.findAll { k, _v -> !(k in ch_chunk_keys) }
            tuple(groupKey(key, meta.hmm_chunks * (meta.read_chunks ?: 1)), domtbl, parquet)
        }
        .groupTuple()
        .map { meta, domtbls, parquets ->
            tuple(meta, domtbls.sort { it.name }, parquets.sort { it.name })
        }

    COMBINEHMMSEARCHTBL(ch_grouped.map { meta, domtbls, _parquets -> tuple(meta, domtbls) })
    MERGEPARQUET(ch_grouped.map { meta, _domtbls, parquets -> tuple(meta, parquets) })
    ch_versions = ch_versions.mix(MERGEPARQUET.out.versions)

    emit:
    domtbl = COMBINEHMMSEARCHTBL.out.concatenated_result // channel: [ val(meta), path(domtbl.gz) ]
    alignments = MERGEPARQUET.out.parquet // channel: [ val(meta), path(parquet) ]
    versions = ch_versions // channel: [ versions.yml ]
}
