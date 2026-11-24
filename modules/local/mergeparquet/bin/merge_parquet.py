#!/usr/bin/env python3
"""Alignment Parquets of one sample into one, sorted by family then read.

A few million domains are read at a time, a range of families from every input, so memory hardly
grows with the sample: 4.6 GB at 24 M domains, where all at once took up to 25 GB on Codon.
"""

import argparse
from collections import Counter

import pyarrow.parquet as pq

ROWS = 2_000_000  # domains per pass, unless one family has more


def family_ranges(files, rows):
    """(first, last) family of each pass."""
    counts = Counter()
    for f in files:
        for batch in pq.ParquetFile(f).iter_batches(columns=["pfam"]):
            for x in batch.column(0).value_counts().to_pylist():
                counts[x["values"]] += x["counts"]
    first, total = None, 0
    for family in sorted(counts):
        if total and total + counts[family] > rows:
            yield first, last
            total = 0
        if not total:
            first = family
        last, total = family, total + counts[family]
    if total:
        yield first, last


def merge(files, out, rows=ROWS):
    schema = pq.read_schema(files[0])
    with pq.ParquetWriter(out, schema, compression="zstd", compression_level=9, use_dictionary=["pfam"]) as writer:
        ranges = list(family_ranges(files, rows))
        for first, last in ranges:
            table = pq.read_table(files, filters=[("pfam", ">=", first), ("pfam", "<=", last)], use_threads=False)
            writer.write_table(table.sort_by([("pfam", "ascending"), ("read", "ascending")]))
        if not ranges:
            writer.write_table(schema.empty_table())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    merge(sorted(args.files), args.out)
