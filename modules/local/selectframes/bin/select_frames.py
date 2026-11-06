#!/usr/bin/env python3
"""Keep ORFs in FragGeneScanRs's frame, plus other ORFs with a hexamer coding score at or above a threshold.

Hexamer log-odds compare ORFs in called frames with ORFs in other frames, from evenly spaced called reads.
Two tables are trained, and each read is scored by the one it did not train. Scores are sums, so they
scale with the training size: the threshold is calibrated for TRAIN_READS, which is deliberately fixed.
"""

import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from math import log

TRAIN_READS = 25000  # per table
ORF_NAME = re.compile(r"^(.*)_frame=(-?\d+)_begin=(\d+)_end=(\d+)$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def open_text(path, mode="rt"):
    return gzip.open(path, mode, compresslevel=1) if str(path).endswith(".gz") else open(path, mode)


def fasta(path):
    """(header line without >, sequence), with a sequence on any number of lines."""
    header, seq = None, []
    with open_text(path) as f:
        for line in f:
            if line[0] == ">":
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def read_calls(path):
    """Called frames per read as (on minus strand, phase), at both ends since a frameshift moves the frame."""
    calls, read = defaultdict(set), None
    with open_text(path) as f:
        for line in f:
            if line[0] == ">":
                read = line[1:].split()[0]
            else:
                start, end, strand = line.split()[:3]
                minus = strand == "-"
                calls[read] |= {(minus, (int(start) - 1) % 3), (minus, int(end) % 3)}
    return calls


def nucleotides(read_seq, header):
    """The read's name, the ORF's frame as in read_calls, and its nucleotides in coding direction."""
    read, frame, begin, end = ORF_NAME.match(header.split()[0]).groups()
    minus, start = int(frame) < 0, int(begin) - 1
    nt = read_seq[start : int(end)]
    return read, (minus, start % 3), nt.translate(COMPLEMENT)[::-1] if minus else nt


def by_read(reads_path, orfs_path):
    """(read name, read sequence, [(ORF header, ORF sequence)]); seqkit keeps ORFs in read order."""
    orfs = fasta(orfs_path)
    pending = next(orfs, None)
    for header, seq in fasta(reads_path):
        read, seq, rows = header.split()[0], seq.upper(), []
        while pending and ORF_NAME.match(pending[0].split()[0])[1] == read:
            rows.append(pending)
            pending = next(orfs, None)
        yield read, seq, rows
    if pending:
        raise ValueError(f"ORF {pending[0]} isn't in read order, or its read is missing")


def hexamers(nt):
    return (nt[i : i + 6] for i in range(0, len(nt) - 5, 3))


def log_odds(coding, other):
    nc, no = sum(coding.values()) + 4096, sum(other.values()) + 4096
    return {h: log((coding[h] + 1) / nc) - log((other[h] + 1) / no) for h in coding.keys() | other.keys()}


def train(reads_path, orfs_path, calls):
    step = max(1, len(calls) // (2 * TRAIN_READS))
    coding, other = [Counter(), Counter()], [Counter(), Counter()]
    called = 0
    for read, seq, rows in by_read(reads_path, orfs_path):
        if read in calls:
            if called % step == 0:
                fold = (called // step) % 2
                for header, _ in rows:
                    _, frame, nt = nucleotides(seq, header)
                    (coding if frame in calls[read] else other)[fold].update(hexamers(nt))
            called += 1
    return [log_odds(coding[1 - i], other[1 - i]) for i in (0, 1)], step


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", required=True, help="The FASTA that seqkit translated, may be gzipped")
    parser.add_argument("--orfs", required=True, help="seqkit translate's output, may be gzipped")
    parser.add_argument("--calls", required=True, help="FragGeneScanRs's .out")
    parser.add_argument("--threshold", type=float, default=-5.0)
    parser.add_argument("--min-calls", type=int, default=5000, help="Below this many called reads, keep all ORFs")
    parser.add_argument("--out", required=True, help="Kept ORFs; gzipped if the name ends in .gz")
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    calls = read_calls(args.calls)
    stats = Counter()
    fallback = len(calls) < args.min_calls
    tables, step = (None, 1) if fallback else train(args.reads, args.orfs, calls)
    called = uncalled = 0
    with open_text(args.out, "wt") as out:
        for read, seq, rows in by_read(args.reads, args.orfs):
            if read in calls:
                fold, called = (called // step) % 2, called + 1
            else:
                fold, uncalled = uncalled % 2, uncalled + 1
            for header, protein in rows:
                stats["orfs_in"] += 1
                stats["residues_in"] += len(protein)
                _, frame, nt = nucleotides(seq, header)
                keep = (
                    fallback
                    or frame in calls.get(read, ())
                    or sum(tables[1 - fold].get(h, 0.0) for h in hexamers(nt)) >= args.threshold
                )
                if keep:
                    out.write(f">{header}\n{protein}\n")
                    stats["orfs_out"] += 1
                    stats["residues_out"] += len(protein)
    stats.update(reads_called=called, reads_uncalled=uncalled)
    with open(args.stats, "w") as f:
        json.dump({**stats, "fallback_kept_all_orfs": fallback, "threshold": args.threshold}, f, indent=1)


if __name__ == "__main__":
    main()
