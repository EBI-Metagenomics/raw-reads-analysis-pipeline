#!/usr/bin/env python3
"""Alignment Parquets of one sample into one, sorted by family then read.

Each input is read once, a row group at a time, and spilled to temporary files that hold a range of
families each, a few million domains. They are then sorted one by one, so memory doesn't grow with
the sample unless one family is larger than a range.
"""

import argparse
import tempfile
from bisect import bisect_left
from collections import Counter
from pathlib import Path

import pyarrow as pa
import pyarrow.ipc
import pyarrow.parquet as pq

ROWS = 2_000_000  # domains sorted at a time, unless one family has more
KEYS = [("pfam", "ascending"), ("read", "ascending")]


def row_groups(files, columns=None):
    # pre_buffer would keep what it has read of the file
    for f in files:
        with pq.ParquetFile(f, pre_buffer=False) as p:
            for i in range(p.num_row_groups):
                yield p.read_row_group(i, columns, use_threads=False)


def range_starts(files, rows):
    """First family of each range."""
    counts = Counter()
    for table in row_groups(files, ["pfam"]):
        for x in table["pfam"].value_counts().to_pylist():
            counts[x["values"]] += x["counts"]
    starts, total = [], 0
    for family in sorted(counts):
        if not starts or total + counts[family] > rows:
            starts.append(family)
            total = 0
        total += counts[family]
    return starts


def spill(files, starts, schema, paths):
    """Every row to the file of its family's range."""
    options = pa.ipc.IpcWriteOptions(compression="zstd", use_threads=False)
    writers = [pa.ipc.new_stream(str(p), schema, options=options) for p in paths]
    for table in row_groups(files):
        table = table.sort_by(KEYS)
        families = table["pfam"]
        edges = [bisect_left(range(len(table)), s, key=lambda i: families[i].as_py()) for s in starts] + [len(table)]
        for writer, start, end in zip(writers, edges, edges[1:]):
            if end > start:
                writer.write_table(table.slice(start, end - start))
    for writer in writers:
        writer.close()


def merge(files, out, rows=ROWS):
    schema = pq.read_schema(files[0])
    starts = range_starts(files, rows)
    with (
        tempfile.TemporaryDirectory(dir=Path(out).parent) as d,
        pq.ParquetWriter(out, schema, compression="zstd", compression_level=9, use_dictionary=[c for c in ("pfam", "pfam_version", "above_ga") if c in schema.names]) as writer,
    ):
        paths = [Path(d) / f"{i}.arrow" for i in range(len(starts))]
        spill(files, starts, schema, paths)
        for path in paths:
            with pa.OSFile(str(path)) as f:
                table = pa.ipc.open_stream(f).read_all()
            path.unlink()
            writer.write_table(table.sort_by(KEYS))
        if not paths:
            writer.write_table(schema.empty_table())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    merge(sorted(args.files), args.out)
