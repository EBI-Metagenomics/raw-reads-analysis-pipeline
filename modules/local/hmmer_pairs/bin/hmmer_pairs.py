#!/usr/bin/env python3
"""HMMER on candidate (model, ORF) pairs only, at Pfam's gathering thresholds.

Writes a domtbl as `hmmsearch --cut_ga --domtblout` does, and the alignments as Parquet if asked to. With
--weak-evalue the Parquet also has the domains under the gathering thresholds; the domtbl never does.

The pairs come from a file of seeds, or with --stream on stdin from kfpfam's `filter join`, run with
FILTER_BINARY=1 beside this script on the same ORFs. A stream is searched as it arrives, a few blocks
of ORFs at a time, so neither the pairs nor the sample's ORFs are ever held or written whole. HMMER
needs little memory itself: what grows is the ORFs, about 1 KB each as Easel sequences with a name to
find them by, for 30 residues.
"""

import argparse
import gzip
import io
import re
import resource
import struct
import sys
import threading
from array import array
from bisect import bisect_right
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager

import pyarrow as pa
import pyarrow.parquet as pq
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile
from pyhmmer.plan7 import Background, HMMFile, Pipeline, Profile

ALPHABET = Alphabet.amino()
ORF_NAME = re.compile(r"^(.*?)_frame=(-?\d+)_begin=(\d+)_end=(\d+)$")
BATCH = 2_000_000  # streamed pairs searched together: at most that many ORFs held, about 1.6 GB
JOB = 512  # pairs per thread job, so that small models don't cost more in hand-over than in HMMER
ROW_GROUP = 65_536  # domains held before they are written
LOCAL = threading.local()
SCHEMA = pa.schema(
    [
        ("read", pa.string()),
        ("frame", pa.int8()),
        ("orf_begin", pa.int32()),
        ("orf_end", pa.int32()),
        ("read_from", pa.int32()),
        ("read_to", pa.int32()),
        ("pfam", pa.string()),
        ("pfam_version", pa.int16()),
        ("above_ga", pa.bool_()),
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


def read_span(frame, begin, end, ali_from, ali_to):
    """The aligned residues' nucleotides on the read: 1-based, inclusive, forward numbering.

    Carried in the output so that anchoring a read to anything else, a protein database say, needs
    no frame arithmetic. `begin` and `end` are the translated stretch as seqkit names it, already
    offset for the frame, and its length is a whole number of codons. A forward frame reads up
    from `begin`, a reverse frame down from `end`, three nucleotides per residue either way.

    Not the mapping in hmmer_domtbl_parse_coverage.py: that one adds the frame offset a second
    time on top of `begin`, and places a reverse frame one nucleotide low. Those shifts cancel
    when it compares hits on one read, which is all it does with them.
    """
    if frame > 0:
        return begin + 3 * (ali_from - 1), begin + 3 * (ali_to - 1) + 2
    return end - 3 * (ali_to - 1) - 2, end - 3 * (ali_from - 1)


def read_seeds(path):
    """Candidates per model name as integer arrays, and the integer of each ORF name."""
    pairs, ids = defaultdict(lambda: array("I")), {}
    with open(path) as f:
        for line in f:
            if line[0] != "#":
                model, target = line.split(None, 2)[:2]
                pairs[model].append(ids.setdefault(target, len(ids)))
    return pairs, ids


@contextmanager
def open_orfs(path):
    with (gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")) as handle:
        with SequenceFile(handle, format="fasta", digital=True, alphabet=ALPHABET) as f:
            yield f


def read_orfs(path, ids):
    """The ORFs named in ids, indexed by their integer."""
    seqs = [None] * len(ids)
    if not ids:
        return seqs
    with open_orfs(path) as f:
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


def domain_rows(hit, domains, pfam, version, above_ga, first=1):
    """A Parquet row for each of a hit's domains, numbered from `first`."""
    read, frame, begin, end = ORF_NAME.match(text(hit.name)).groups()
    frame, begin, end = int(frame), int(begin), int(end)
    for n, domain in enumerate(domains, first):
        ali = domain.alignment
        read_from, read_to = read_span(frame, begin, end, ali.target_from, ali.target_to)
        yield (
            read, frame, begin, end, read_from, read_to, pfam, version, above_ga, n,
            hit.score, domain.score, hit.evalue, domain.i_evalue,
            ali.hmm_from, ali.hmm_to, ali.target_from, ali.target_to,
            domain.env_from, domain.env_to,
            cigar(ali.hmm_sequence, ali.target_sequence),
            ali.target_sequence.replace("-", "").upper(),
            ali.posterior_probabilities,
        )  # fmt: skip


def pipeline(name, **options):
    """This thread's pipeline of that name, cleared. One is kept a thread, as building one is slow."""
    if not hasattr(LOCAL, name):
        setattr(LOCAL, name, Pipeline(ALPHABET, **options))
    p = getattr(LOCAL, name)
    p.clear()
    return p


def search(query, block, num_models, weak_evalue=0):
    """One model on its candidates: domtbl bytes, and one row per domain.

    At gathering thresholds, and with `weak_evalue` also the domains under them whose ORF and own E-values are
    at most that. Every candidate is then searched without a cut first, and at the gathering thresholds only the ORFs that
    search reports: it stops an ORF at the same filters, so the domtbl is as it is without `weak_evalue`.
    """
    accession, _, version = text(query.accession).partition(".")
    version = int(version or 0)
    weak = []
    if weak_evalue > 0:
        every = pipeline("every", Z=num_models, E=float("inf"), domE=float("inf")).search_hmm(query, block)
        found = {text(hit.name) for hit in every}
        block = DigitalSequenceBlock(ALPHABET, [seq for seq in block if text(seq.name) in found])
        weak = [hit for hit in every if hit.evalue <= weak_evalue]
    hits = pipeline("gathering", Z=num_models, bit_cutoffs="gathering").search_hmm(query, block)
    table = io.BytesIO()
    hits.write(table, format="domains", header=False)
    rows, above = [], {}
    for hit in hits:
        if hit.reported:
            domains = [d for d in hit.domains if d.reported]
            above[text(hit.name)] = {(d.env_from, d.env_to) for d in domains}
            rows.extend(domain_rows(hit, domains, accession, version, True))
    for hit in weak:
        done = above.get(text(hit.name), set())
        domains = [d for d in hit.domains if d.i_evalue <= weak_evalue and (d.env_from, d.env_to) not in done]
        rows.extend(domain_rows(hit, domains, accession, version, False, len(done) + 1))
    return table.getvalue(), rows


class Output:
    """The domtbl, gzipped, and the alignments, zstd Parquet (none if its path is None), written as they come.

    A domain takes about 700 bytes as a Python tuple and 55 in the file, so rows are held only until
    they make a row group. Each group is sorted by family then read; MERGEPARQUET sorts the whole.
    """

    def __init__(self, domtbl, parquet):
        self.domtbl, self.rows = gzip.open(domtbl, "wb"), []
        self.parquet = parquet and pq.ParquetWriter(
            parquet, SCHEMA, compression="zstd", compression_level=9, use_dictionary=["pfam", "pfam_version", "above_ga"]
        )

    def add(self, table, rows):
        self.domtbl.write(table)
        if self.parquet:
            self.rows.extend(rows)
            if len(self.rows) >= ROW_GROUP:
                self.flush()

    def flush(self):
        columns = [list(c) for c in zip(*self.rows)] or [[] for _ in SCHEMA.names]
        table = pa.Table.from_arrays([pa.array(c, f.type) for c, f in zip(columns, SCHEMA)], schema=SCHEMA)
        self.parquet.write_table(table.sort_by([("pfam", "ascending"), ("read", "ascending")]))
        self.rows = []

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.domtbl.close()
        if self.parquet:
            self.flush()
            self.parquet.close()


def run(jobs, num_models, threads, pool, output, weak_evalue=0):
    """Every (model, its candidate ORFs) of jobs, written in their order so output doesn't depend on threads.

    A model may be in jobs once: an optimized profile is reconfigured for each target, so two threads
    can't share one.
    """

    def work(group):
        return [search(query, block, num_models, weak_evalue) for query, block in group]

    def collect(future):
        for table, rows in future.result():
            output.add(table, rows)

    pending, group, n = deque(), [], 0
    for query, seqs in jobs:
        # Built here, not in the worker: blocks built across 6 threads ran about 40x slower
        group.append((query, DigitalSequenceBlock(ALPHABET, seqs)))
        n += len(seqs)
        if n >= JOB:
            pending.append(pool.submit(work, group))
            group, n = [], 0
            if len(pending) >= 4 * threads:
                collect(pending.popleft())
    if group:
        pending.append(pool.submit(work, group))
    while pending:
        collect(pending.popleft())


def seeded(pairs, seqs, hmm_path, whole=False):
    """(model, its candidate ORFs) from seeds, in HMM file order.

    Seeds name models, and a chunk of a chunked database skips the names of the other chunks. When
    the file is the whole database (`whole`), a name it lacks means the seeds came from another Pfam release:
    skipping those would lose every renamed family without a word, so it is an error.
    """
    names = set()
    with HMMFile(hmm_path) as hmms:
        for hmm in hmms:
            name = text(hmm.name)
            names.add(name)
            targets = pairs.get(name)
            if targets:
                yield hmm, [seqs[i] for i in sorted(set(targets))]
    unknown = sorted(pairs.keys() - names) if whole else []
    if unknown:
        raise ValueError(
            f"{len(unknown)} seeded models are not in {hmm_path}, {', '.join(unknown[:5])} among them: "
            "the seeds and the HMM file are from different Pfam releases"
        )


def profiles(hmm_paths, models_tsv):
    """The optimized profile of each model of kfpfam's models.tsv, in its order, which is what a stream's
    model numbers index. Held for the run, as a stream comes back to every model batch after batch."""
    background, found = Background(ALPHABET), {}
    for path in hmm_paths:
        with HMMFile(path) as hmms:
            for hmm in hmms:
                profile = Profile(hmm.M, ALPHABET)
                profile.configure(hmm, background)
                found[text(hmm.name)] = profile.to_optimized()
    with open(models_tsv) as f:
        names = [line.rstrip("\n").split("\t")[1] for line in f][1:]
    if set(names) != found.keys():
        odd = sorted(set(names) ^ found.keys())
        raise ValueError(
            f"{len(odd)} models are in {models_tsv} or the HMM files but not both, {', '.join(odd[:5])} among them: "
            "they are from different Pfam releases"
        )
    return [found[n] for n in names]


def read_blocks(stream, orfs):
    """`filter join`'s blocks: their ORFs, read forward from the FASTA it read, and their pairs, as one flat
    array of ORF in the block, model, ORF, model... in model order."""

    def take(size):
        data = stream.read(size)
        if len(data) < size:
            raise ValueError("the filter's output stops inside a block")
        return data

    while head := stream.read(24):
        n, residues, size = struct.unpack("<3Q", head + take(24 - len(head)))
        pairs = array("I")
        pairs.frombytes(take(size))
        block = orfs.read_block(sequences=n)
        if (len(block), sum(map(len, block))) != (n, residues):
            raise ValueError("the filter's blocks are not of these ORFs")
        yield block, pairs
    if orfs.read() is not None:
        raise ValueError("the filter's output stops before the last ORF")


def batches(blocks, size):
    """Candidate ORFs by model, of as many blocks as make `size` pairs. Only ORFs with a candidate are kept."""
    targets, n = defaultdict(list), 0
    for block, pairs in blocks:
        orfs, models = pairs[::2], pairs[1::2]
        i = 0
        while i < len(models):
            end = bisect_right(models, models[i], i)
            targets[models[i]].extend([block[k] for k in orfs[i:end]])
            i = end
        n += len(models)
        if n >= size:
            yield targets
            targets, n = defaultdict(list), 0
    if targets:
        yield targets


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--orfs", required=True, help="ORF FASTA, may be gzipped")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--seeds", help="nail's --seeds-out: model name, then ORF name")
    source.add_argument("--stream", metavar="MODELS_TSV", help="pairs on stdin, of the models in kfpfam's models.tsv")
    parser.add_argument("--hmm", required=True, nargs="+", help="one file, or with --stream every chunk")
    parser.add_argument("--num-models", type=int, required=True, help="Z, as in the pipeline's hmmsearch")
    parser.add_argument("--whole", action="store_true", help="--hmm is the whole database, not a chunk of it")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--domtbl", required=True, help="Gzipped output")
    parser.add_argument("--parquet", help="Alignments, not written if this is left out")
    parser.add_argument("--weak-evalue", type=float, default=0, metavar="E",
                        help="domains under the gathering thresholds too, in the Parquet only, at this E-value or less (0: none)")
    args = parser.parse_args()
    if args.weak_evalue < 0:
        parser.error("--weak-evalue is an E-value, or 0 for none")
    weak = args.weak_evalue if args.parquet else 0
    if args.seeds and len(args.hmm) > 1:
        parser.error("--seeds takes one HMM file")
    with Output(args.domtbl, args.parquet) as output, ThreadPoolExecutor(args.threads) as pool:
        if args.seeds:
            pairs, ids = read_seeds(args.seeds)
            seqs = read_orfs(args.orfs, ids)
            del ids
            run(seeded(pairs, seqs, args.hmm[0], args.whole), args.num_models, args.threads, pool, output, weak)
        else:
            optimized = profiles(args.hmm, args.stream)
            with open_orfs(args.orfs) as orfs:
                for targets in batches(read_blocks(sys.stdin.buffer, orfs), BATCH):
                    jobs = ((optimized[m], targets[m]) for m in sorted(targets))
                    run(jobs, args.num_models, args.threads, pool, output, weak)
    used = resource.getrusage(resource.RUSAGE_SELF)
    print(f"HMMER: {used.ru_utime + used.ru_stime:.0f} core-s", file=sys.stderr)


if __name__ == "__main__":
    main()
