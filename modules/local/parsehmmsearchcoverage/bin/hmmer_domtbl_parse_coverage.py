"""HMM model coverages of metagenomic reads, from hmmsearch's domain table.

A read's overlapping hits are resolved greedily, best E-value first, so all of its hits are needed at once and the
table isn't sorted by read. Hits are split by read into temporary files, which are then resolved one at a time: memory
is that of one file's hits (1/64 of them), where holding the whole table took about 1 GB per million reads.
"""

import argparse
import gc
import json
import marshal
import sys
import tempfile
from collections import defaultdict
from itertools import accumulate

PARTS = 64
BUFFER = 4096  # hits of a part held before they are written


def parse(lines):
    """A tuple per hit: read, line number, E-value, accession, model length, start and end on the read, and on the model.

    The read's coordinates are only compared with each other, and are 0-based with the end excluded.
    """
    for n, line in enumerate(lines):
        if line.startswith(b"#"):
            continue
        f = line.split(None, 19)
        read, _, orf = f[0].rpartition(b"_frame=")
        frame, begin, end = orf.split(b"_")
        frame, ali_from, ali_to = int(frame), int(f[17]), int(f[18]) + 1
        if frame < 0:
            offset = int(end[4:]) + 1 - frame
            start, stop = offset - 3 * ali_to, offset - 3 * ali_from
        else:
            offset = int(begin[6:]) - 5 + frame
            start, stop = offset + 3 * ali_from, offset + 3 * ali_to
        yield read, n, float(f[6]), f[4], int(f[5]), start, stop, int(f[15]), int(f[16])


def split(hits, files):
    buffers = [[] for _ in files]
    for hit in hits:
        part = hash(hit[0]) % len(files)
        buffers[part].append(hit)
        if len(buffers[part]) == BUFFER:
            marshal.dump(buffers[part], files[part])
            buffers[part] = []
    for buffer, file in zip(buffers, files):
        marshal.dump(buffer, file)
        marshal.dump([], file)  # the end
        file.seek(0)


def top_hits(hits):
    """A read's hits that don't overlap a better one on it, best first."""
    if len(hits) == 1:
        return hits
    top = []
    for hit in sorted(hits, key=lambda h: h[2]):
        if all(hit[6] <= t[5] or t[6] <= hit[5] for t in top):
            top.append(hit)
    return top


def count(files):
    """Reads mapped, and per model: reads, depth at each position, and its place in the order models are first seen."""
    reads_mapped, reads, steps, first = 0, defaultdict(int), {}, {}
    for file in files:
        by_read = defaultdict(list)
        while chunk := marshal.load(file):
            for hit in chunk:
                by_read[hit[0]].append(hit)
        file.close()
        reads_mapped += len(by_read)
        for hits in by_read.values():
            line, seen = hits[0][1], set()
            for rank, hit in enumerate(top_hits(hits)):
                _, _, _, model, length, _, _, hmm_from, hmm_to = hit
                if model not in steps:
                    steps[model] = [0] * (length + 1)
                    first[model] = line, rank
                elif (line, rank) < first[model]:
                    first[model] = line, rank
                # depth changes here: positions are 1-based with the end included
                steps[model][hmm_from - 1] += 1
                steps[model][hmm_to] -= 1
                if model not in seen:
                    seen.add(model)
                    reads[model] += 1
    depth = {model: list(accumulate(steps[model][:-1])) for model in sorted(steps, key=first.get)}
    return reads_mapped, reads, depth


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_fp", default="-", help="Domain table, '-' for stdin (default)")
    parser.add_argument("-o", "--output_fp", default="-", help="Output TSV, '-' for stdout (default)")
    parser.add_argument("-s", "--stats_output_fp", default="", help="Output JSON of mapping stats")
    args = parser.parse_args()
    gc.disable()  # nothing is cyclic, and it walks the hits held

    files = [tempfile.TemporaryFile(dir=".") for _ in range(PARTS)]
    with sys.stdin.buffer if args.input_fp == "-" else open(args.input_fp, "rb") as lines:
        split(parse(lines), files)
    reads_mapped, reads, coverage = count(files)

    rows = [(model, sum(d) / len(d), sum(v > 0 for v in d) / len(d)) for model, d in coverage.items() if d]
    with sys.stdout if args.output_fp == "-" else open(args.output_fp, "wt") as out:
        out.write("# function\tread_count\tcoverage_depth\tcoverage_breadth\n")
        for model, depth, breadth in sorted(rows, key=lambda r: -r[1]):
            out.write(f"{model.decode()}\t{reads[model]}\t{depth}\t{breadth}\n")

    if args.stats_output_fp:
        stats = {"reads_mapped": reads_mapped, "hmm_count": len(reads), "read_hit_count": sum(reads.values())}
        with open(args.stats_output_fp, "wt") as f:
            json.dump(stats, f)


if __name__ == "__main__":
    main()
