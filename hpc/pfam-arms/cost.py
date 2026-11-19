#!/usr/bin/env python3
"""Core-hours and peak memory of each arm's Pfam steps, from Nextflow traces.

Usage: cost.py ARM=TRACE.txt ...
Core-hours are realtime x allocated CPUs, for tasks that weren't cached in that arm's run.
"""

import csv
import re
import sys
from collections import defaultdict

STEP = re.compile(r"FRAGGENESCANRS|SELECTFRAMES|NAIL_SEED|HMMER_PAIRS|HMMER_HMMSEARCH|SEQKIT_TRANSLATE|PARSEHMMSEARCHCOVERAGE")
UNITS = {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}
BYTES = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}


def seconds(text):
    return sum(float(n) * UNITS[u] for n, u in re.findall(r"([\d.]+)\s*(ms|s|m|h|d)", text))


def gigabytes(text):
    m = re.match(r"([\d.]+)\s*([KMGT]?B)", text)
    return float(m[1]) * BYTES[m[2]] / 1e9 if m else 0.0


def main():
    print("arm\tsample\tstep\ttasks\tcore_hours\tpeak_rss_gb")
    for spec in sys.argv[1:]:
        arm, path = spec.split("=", 1)
        rows = defaultdict(lambda: [0, 0.0, 0.0])
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                step = STEP.search(r["process"])
                if step and r["status"] == "COMPLETED":
                    row = rows[(r["tag"], step[0])]
                    row[0] += 1
                    row[1] += seconds(r["realtime"]) * int(r["cpus"]) / 3600
                    row[2] = max(row[2], gigabytes(r["peak_rss"]))
        for (sample, step), (tasks, hours, rss) in sorted(rows.items()):
            print(f"{arm}\t{sample}\t{step}\t{tasks}\t{hours:.3f}\t{rss:.2f}")


if __name__ == "__main__":
    main()
