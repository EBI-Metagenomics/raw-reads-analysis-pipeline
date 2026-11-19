#!/usr/bin/env python3
"""Compare Pfam profiles from different arms of one sample, each against the first.

Usage: compare_profiles.py SAMPLE BASE_NAME=BASE.tsv[.gz] NAME=ARM.tsv[.gz] ... [--families OUT.tsv]
Prints one row per arm: reads and families, read-count differences, and the worst families among
those with at least 10 reads in the base. Coverage breadth is compared on families in both.
"""

import gzip
import sys


def profile(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        rows = (x.rstrip("\n").split("\t") for x in f if x[0] != "#")
        return {c[0].split(".")[0]: (int(c[1]), float(c[3])) for c in rows}


def compare(sample, base, arm):
    fams = base.keys() | arm.keys()
    moved = sum(abs(base.get(k, (0, 0))[0] - arm.get(k, (0, 0))[0]) for k in fams)
    total = sum(v[0] for v in base.values())
    lost = [k for k in base if k not in arm]
    common = [k for k in fams if k in base and k in arm]
    breadth = sum(abs(base[k][1] - arm[k][1]) for k in common) / max(1, len(common))
    big = [k for k in base if base[k][0] >= 10]
    ratios = sorted((arm.get(k, (0, 0))[0] / base[k][0], k, base[k][0]) for k in big)
    below = sum(r < 0.8 for r, _, _ in ratios)
    return {
        "reads": sum(v[0] for v in arm.values()),
        "families": len(arm),
        "abs_read_diff_pct": round(100 * moved / max(1, total), 2),
        "families_lost": len(lost),
        "reads_in_lost_families": sum(base[k][0] for k in lost),
        "families_gained": sum(k not in base for k in arm),
        "mean_breadth_diff": round(breadth, 4),
        "families_10plus": len(big),
        "of_those_below_0.8": below,
        "worst": ";".join(f"{k}:{r:.2f}({n})" for r, k, n in ratios[:3]),
    }


def main():
    args = sys.argv[1:]
    out = None
    if "--families" in args:
        i = args.index("--families")
        out = args[i + 1]
        del args[i : i + 2]
    sample, specs = args[0], [a.split("=", 1) for a in args[1:]]
    profiles = {name: profile(path) for name, path in specs}
    base_name = specs[0][0]
    base = profiles[base_name]
    keys = None
    for name, _ in specs[1:]:
        row = compare(sample, base, profiles[name])
        if keys is None:
            keys = list(row)
            print("sample\tarm\t" + "\t".join(keys))
        print(f"{sample}\t{name}\t" + "\t".join(str(row[k]) for k in keys))
    print(f"{sample}\t{base_name} (base)\treads {sum(v[0] for v in base.values())}\tfamilies {len(base)}")
    if out:
        fams = sorted(set().union(*profiles.values()))
        with open(out, "w") as f:
            f.write("sample\tfamily\t" + "\t".join(n for n, _ in specs) + "\n")
            for k in fams:
                f.write(f"{sample}\t{k}\t" + "\t".join(str(profiles[n].get(k, (0, 0))[0]) for n, _ in specs) + "\n")


if __name__ == "__main__":
    main()
