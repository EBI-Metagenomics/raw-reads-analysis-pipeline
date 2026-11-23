#!/usr/bin/env python3
"""Core-hours and peak memory of each arm's Pfam steps, from Nextflow traces.

Usage: cost.py ARM=TRACE.txt ...
cpu_hours is realtime x %cpu, what the tasks used; core_hours is realtime x allocated CPUs, empty
when the trace has no cpus column. Only tasks that weren't cached in that arm's run are counted. `failed` is
the attempts that failed: a killed one leaves no CPU time in the trace, so their cost is not here.
"""

import csv
import re
import sys
from collections import defaultdict

STEP = re.compile(
    r"FRAGGENESCANRS|SELECTFRAMES|NAIL_SEED|KFPFAM_FILTER|HMMER_PAIRS|HMMER_HMMSEARCH|SEQKIT_TRANSLATE|"
    r"COMBINEHMMSEARCHTBL|MERGEPARQUET|PARSEHMMSEARCHCOVERAGE"
)
UNITS = {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}
BYTES = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}


def seconds(text):
    return sum(float(n) * UNITS[u] for n, u in re.findall(r"([\d.]+)\s*(ms|s|m|h|d)", text))


def gigabytes(text):
    m = re.match(r"([\d.]+)\s*([KMGT]?B)", text)
    return float(m[1]) * BYTES[m[2]] / 1e9 if m else 0.0


def main():
    print("arm\tsample\tstep\ttasks\tcpu_hours\tcore_hours\tpeak_rss_gb\tfailed")
    for spec in sys.argv[1:]:
        arm, path = spec.split("=", 1)
        rows, has_cpus = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0]), False
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                # Nextflow's default fields have neither process nor tag: both are in name, "A:B:STEP (tag)"
                process, _, tag = r["name"].partition(" (")
                step = STEP.search(process.rsplit(":", 1)[-1])
                if not step or r["status"] not in ("COMPLETED", "FAILED"):
                    continue
                row = rows[(tag.rstrip(")"), step[0])]
                hours = seconds(r["realtime"]) / 3600
                cpu = r["%cpu"].rstrip("%")
                used = hours * float(cpu) / 100 if cpu not in ("", "-") else 0.0
                if r["status"] == "FAILED":
                    row[4] += 1
                    continue
                row[0] += 1
                row[1] += used
                row[2] += hours * int(r["cpus"]) if r.get("cpus") else 0.0
                row[3] = max(row[3], gigabytes(r["peak_rss"]))
                has_cpus = has_cpus or bool(r.get("cpus"))
        for (sample, step), (tasks, used, core, rss, failed) in sorted(rows.items()):
            core = f"{core:.3f}" if has_cpus else ""
            print(f"{arm}\t{sample}\t{step}\t{tasks}\t{used:.3f}\t{core}\t{rss:.2f}\t{failed}")


if __name__ == "__main__":
    main()
