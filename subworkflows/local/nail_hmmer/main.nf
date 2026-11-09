include { NAIL_SEED } from '../../../modules/local/nail_seed/main'
include { HMMER_PAIRS } from '../../../modules/local/hmmer_pairs/main'

workflow NAIL_HMMER {
    take:
    orfs // channel: [ val(meta), path(fasta) ]
    pfam_db // channel: path(hmm)

    main:
    // One whole Pfam-A.hmm. PROFILE_HMMSEARCH_PFAM rejects a chunked database before reaching
    // here, so this channel has exactly one element; .first() makes it a value channel that
    // every sample can consume.
    ch_hmm = pfam_db.first()

    // nail stops on an empty target file; hmmsearch's chunking also skips samples with no ORFs
    ch_orfs = orfs.filter { _meta, fasta -> fasta.countFasta(limit: 1) > 0 }

    NAIL_SEED(ch_orfs, ch_hmm)

    HMMER_PAIRS(
        ch_orfs.join(NAIL_SEED.out.seeds),
        ch_hmm,
        params.databases.pfam.variables.num_models,
    )

    emit:
    domtbl = HMMER_PAIRS.out.domtbl // channel: [ val(meta), path(domtbl.gz) ]
    alignments = HMMER_PAIRS.out.alignments // channel: [ val(meta), path(parquet) ]
    versions = NAIL_SEED.out.versions.mix(HMMER_PAIRS.out.versions) // channel: [ versions.yml ]
}
