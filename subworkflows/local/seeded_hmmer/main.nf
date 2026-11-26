include { NAIL_SEED } from '../../../modules/local/nail_seed/main'
include { KFPFAM_HMMER } from '../../../modules/local/kfpfam_hmmer/main'
include { HMMER_PAIRS } from '../../../modules/local/hmmer_pairs/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'
include { MERGEPARQUET } from '../../../modules/local/mergeparquet/main'

// Back to one domtbl or one Parquet per sample: a task per (read chunk, Pfam chunk), so the group is
// as big as their product. A single chunk of each is grouped too, so the chunked and unchunked paths
// stay the same shape. Sorted so the concatenated domtbl does not depend on the order chunks finish in.
def bySample(ch) {
    return ch
        .map { meta, output ->
            def key = meta.findAll { k, _v -> !(k in ['hmm_chunk', 'hmm_chunks', 'read_chunk', 'read_chunks']) }
            tuple(groupKey(key, (meta.hmm_chunks ?: 1) * (meta.read_chunks ?: 1)), output)
        }
        .groupTuple()
        .map { meta, outputs -> tuple(meta, outputs.sort { it.name }) }
}

workflow SEEDED_HMMER {
    take:
    orfs // channel: [ val(meta), path(fasta) ]
    pfam_db // channel: path(hmm), one per chunk of a chunked Pfam, else a single element

    main:
    // Both engines stop on an empty target file; hmmsearch's chunking also skips samples with no ORFs
    ch_orfs = orfs.filter { _meta, fasta -> fasta.countFasta(limit: 1) > 0 }

    ch_versions = Channel.empty()
    // from the command line it is the text 'false' under Nextflow 26's parser, which is true
    def write_alignments = params.pfam_alignments.toString().toBoolean()
    def weak_evalue = params.pfam_weak_evalue.toString().toDouble()

    if (params.pfam_search_engine == 'kfpfam') {
        // kfpfam's index covers the whole database, so this is one task per read chunk whatever the
        // Pfam chunking is, and it filters and runs HMMER together: the pairs go through a pipe
        KFPFAM_HMMER(
            ch_orfs,
            pfam_db.toSortedList(),
            file(params.kfpfam_prep, checkIfExists: true),
            file(params.kfpfam_index, checkIfExists: true),
            params.kfpfam_bar,
            params.kfpfam_gate,
            params.kfpfam_min_gapped,
            params.kfpfam_block,
            params.databases.pfam.variables.num_models,
            write_alignments,
            weak_evalue,
            params.kfpfam_forward_slack,
        )
        ch_versions = ch_versions.mix(KFPFAM_HMMER.out.versions)
        ch_domtbl = KFPFAM_HMMER.out.domtbl
        ch_alignments = KFPFAM_HMMER.out.alignments
    }
    else {
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
            write_alignments,
            weak_evalue,
        )
        ch_versions = ch_versions.mix(NAIL_SEED.out.versions, HMMER_PAIRS.out.versions)
        ch_domtbl = HMMER_PAIRS.out.domtbl
        ch_alignments = HMMER_PAIRS.out.alignments
    }

    COMBINEHMMSEARCHTBL(bySample(ch_domtbl))
    // No input, so no task, without pfam_alignments
    MERGEPARQUET(bySample(ch_alignments))
    ch_versions = ch_versions.mix(MERGEPARQUET.out.versions)

    emit:
    domtbl = COMBINEHMMSEARCHTBL.out.concatenated_result // channel: [ val(meta), path(domtbl.gz) ]
    alignments = MERGEPARQUET.out.parquet // channel: [ val(meta), path(parquet) ]
    versions = ch_versions // channel: [ versions.yml ]
}
