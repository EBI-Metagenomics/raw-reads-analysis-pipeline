#!/usr/bin/env python3
"""HMMER on candidate (model, ORF) pairs only, at Pfam's gathering thresholds.

Writes a domtbl as `hmmsearch --cut_ga --domtblout` does, and the alignments as Parquet.
"""

import argparse
import gzip
import io
import re
from array import array
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor

import pyarrow as pa
import pyarrow.parquet as pq
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile
from pyhmmer.plan7 import HMMFile, Pipeline

ALPHABET = Alphabet.amino()
ORF_NAME = re.compile(r"^(.*?)_frame=(-?\d+)_begin=(\d+)_end=(\d+)$")
SCHEMA = pa.schema(
    [
        ("read", pa.string()),
        ("frame", pa.int8()),
        ("orf_begin", pa.int32()),
        ("orf_end", pa.int32()),
        ("pfam", pa.string()),
        ("domain", pa.int16()),
        ("seq_score", pa.float32()),
        ("dom_score", pa.float32()),
        ("seq_evalue", pa.float64()),
        ("i_evalue", pa.float64()),
        ("hmm_from", pa.int32()),
        ("hmm_to", pa.int32()),
        ("ali_from", pa.int32()),
        ("ali_to", pa.int32()),
        ("env_from", pa.int32()),
        ("env_to", pa.int32()),
        ("cigar", pa.string()),
        ("residues", pa.string()),
        ("pp", pa.string()),
    ]
)


def text(value):
    return value.decode() if isinstance(value, bytes) else value


def read_seeds(path):
    """Candidates per model name as integer arrays, and the integer of each ORF name."""
    pairs, ids = defaultdict(lambda: array("I")), {}
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                model, target = line.split(None, 2)[:2]
                pairs[model].append(ids.setdefault(target, len(ids)))
    return pairs, ids


def read_orfs(path, ids):
    """The ORFs named in ids, indexed by their integer."""
    seqs = [None] * len(ids)
    if not ids:
        return seqs
    with (gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")) as handle:
        with SequenceFile(handle, format="fasta", digital=True, alphabet=ALPHABET) as f:
            for seq in f:
                i = ids.get(text(seq.name))
                if i is not None:
                    seqs[i] = seq
    missing = seqs.count(None)
    if missing:
        raise ValueError(f"{missing} candidate ORFs are not in {path}")
    return seqs


def cigar(hmm_line, target_line):
    """M, I and D against the HMM, from HMMER's aligned lines."""
    ops = "".join("I" if h == "." else "D" if t == "-" else "M" for h, t in zip(hmm_line, target_line))
    return "".join(f"{len(m[0])}{m[1]}" for m in re.finditer(r"(.)\1*", ops))


def search(hmm, block, num_models):
    """One model on its candidates: domtbl bytes, and one row per reported domain."""
    hits = Pipeline(ALPHABET, Z=num_models, bit_cutoffs="gathering").search_hmm(hmm, block)
    table = io.BytesIO()
    hits.write(table, format="domains", header=False)
    pfam = text(hmm.accession).split(".")[0]
    rows = []
    for hit in hits:
        if not hit.reported:
            continue
        name = text(hit.name)
        read, frame, begin, end = ORF_NAME.match(name).groups()
        for n, domain in enumerate((d for d in hit.domains if d.reported), 1):
            ali = domain.alignment
            rows.append(
                (
                    read, int(frame), int(begin), int(end), pfam, n,
                    hit.score, domain.score, hit.evalue, domain.i_evalue,
                    ali.hmm_from, ali.hmm_to, ali.target_from, ali.target_to,
                    domain.env_from, domain.env_to,
                    cigar(ali.hmm_sequence, ali.target_sequence),
                    ali.target_sequence.replace("-", "").upper(),
                    ali.posterior_probabilities,
                )
            )
    return table.getvalue(), rows


def run(pairs, seqs, hmm_path, num_models, threads, domtbl_path, parquet_path):
    """Every model with candidates, written in HMM file order so output doesn't depend on threads."""
    rows = []
    with gzip.open(domtbl_path, "wb") as domtbl, ThreadPoolExecutor(threads) as pool, HMMFile(hmm_path) as hmms:

        def collect(future):
            table, found = future.result()
            domtbl.write(table)
            rows.extend(found)

        pending = deque()
        for hmm in hmms:
            targets = pairs.get(text(hmm.name))
            if targets:
                # Built here, not in the worker: blocks built across 6 threads ran about 40x slower
                block = DigitalSequenceBlock(ALPHABET, [seqs[i] for i in sorted(set(targets))])
                pending.append(pool.submit(search, hmm, block, num_models))
                if len(pending) >= 4 * threads:
                    collect(pending.popleft())
        while pending:
            collect(pending.popleft())
    columns = [list(c) for c in zip(*rows)] or [[] for _ in SCHEMA.names]
    alignments = pa.Table.from_arrays([pa.array(c, f.type) for c, f in zip(columns, SCHEMA)], schema=SCHEMA)
    alignments = alignments.sort_by([("pfam", "ascending"), ("read", "ascending")])
    pq.write_table(alignments, parquet_path, compression="zstd", compression_level=9, use_dictionary=["pfam"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orfs", required=True, help="ORF FASTA, may be gzipped")
    parser.add_argument("--seeds", required=True, help="nail's --seeds-out: model name, then ORF name")
    parser.add_argument("--hmm", required=True)
    parser.add_argument("--num-models", type=int, required=True, help="Z, as in the pipeline's hmmsearch")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--domtbl", required=True, help="Gzipped output")
    parser.add_argument("--parquet", required=True)
    args = parser.parse_args()
    pairs, ids = read_seeds(args.seeds)
    seqs = read_orfs(args.orfs, ids)
    del ids
    run(pairs, seqs, args.hmm, args.num_models, args.threads, args.domtbl, args.parquet)


if __name__ == "__main__":
    main()
