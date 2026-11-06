#!/usr/bin/env python3
"""Development check: nail's alignments against HMMER_PAIRS's, and whether nail alone would do.

GA is called from nail's scores as `--cut_ga` does: the pair's best score at GA1, and the alignment's at GA2.
Usage: compare_alignments.py NAIL.ali PAIRS.parquet PAIRS.domtbl.gz Pfam-A.hmm PARSER.py OUT_DIR
"""

import gzip
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq

LINE = re.compile(r"^\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$")


def cigar(hmm_line, target_line):
    ops = "".join("I" if h in ".-" else "D" if t == "-" else "M" for h, t in zip(hmm_line, target_line))
    return "".join(f"{len(m[0])}{m[1]}" for m in re.finditer(r"(.)\1*", ops))


def columns(hmm_from, ali_from, cig):
    """The (HMM position, ORF position) pairs an alignment aligns."""
    h, t, out = hmm_from, ali_from, set()
    for n, op in re.findall(r"(\d+)([MID])", cig):
        n = int(n)
        if op == "M":
            out.update((h + i, t + i) for i in range(n))
            h, t = h + n, t + n
        elif op == "D":
            h += n
        else:
            t += n
    return out


def read_models(path):
    models, name, acc, length = {}, None, None, None
    with open(path) as f:
        for line in f:
            key = line[:5].strip()
            if key == "ACC":
                acc = line.split()[1]
            elif key == "LENG":
                length = int(line.split()[1])
            elif key == "GA":
                ga = line.split()[1:3]
                models[acc.split(".")[0]] = (float(ga[0]), float(ga[1].rstrip(";")), length, acc)
    return models


def read_nail(path):
    """Alignments per (ORF, Pfam accession)."""
    pairs = defaultdict(list)
    rec = None
    with open(path) as f:
        for line in f:
            if line.startswith("query:"):
                rec = {"q": line.split()[1], "h": "", "t": "", "hf": None, "tf": None}
            elif rec is None:
                continue
            elif line.startswith("target:"):
                rec["orf"] = line.split()[1]
            elif line.startswith("score:"):
                rec["score"] = float(line.split()[1])
            elif line.startswith("E-value:"):
                rec["e"] = line.split()[1]
            elif line.startswith("//"):
                pairs[(rec["orf"], rec["q"].split(".")[0])].append(
                    (rec["score"], rec["e"], rec["hf"], rec["ht"], rec["tf"], rec["tt"], cigar(rec["h"], rec["t"]))
                )
                rec = None
            else:
                m = LINE.match(line)
                if m and m[1] == rec["q"]:
                    rec["h"] += m[3]
                    rec["ht"] = int(m[4])
                    rec["hf"] = rec["hf"] or int(m[2])
                elif m and m[1] == rec["orf"]:
                    rec["t"] += m[3]
                    rec["tt"] = int(m[4])
                    rec["tf"] = rec["tf"] or int(m[2])
    return pairs


def profile(domtbl, parser, stem):
    cat = f"gunzip -c {domtbl}" if str(domtbl).endswith(".gz") else f"cat {domtbl}"
    subprocess.run(f"{cat} | python3 {parser} -i - -o {stem}.tsv -s {stem}.json", shell=True, check=True)
    with open(f"{stem}.tsv") as f:
        return {x.split("\t")[0].split(".")[0]: int(x.split("\t")[1]) for x in f if x[0] != "#"}


def main():
    nail_ali, parquet, domtbl, hmm, parser, out = sys.argv[1:7]
    out = Path(out)
    models = read_models(hmm)
    nail = read_nail(nail_ali)
    print("nail alignments", sum(map(len, nail.values())), "pairs", len(nail))
    called = {}
    for k, alis in nail.items():
        ga = models.get(k[1])
        if ga and max(a[0] for a in alis) >= ga[0]:
            kept = [a for a in alis if a[0] >= ga[1]]
            if kept:
                called[k] = kept
    hmmer = defaultdict(list)
    for r in pq.read_table(parquet).to_pylist():
        k = (f"{r['read']}_frame={r['frame']}_begin={r['orf_begin']}_end={r['orf_end']}", r["pfam"])
        hmmer[k].append(r)

    both, only_n, only_h = called.keys() & hmmer.keys(), called.keys() - hmmer.keys(), hmmer.keys() - called.keys()
    union = len(both) + len(only_n) + len(only_h)
    differ = (len(only_n) + len(only_h)) / union
    print(f"GA pairs: HMMER {len(hmmer)}, nail {len(called)}, both {len(both)}, nail only {len(only_n)}, HMMER only {len(only_h)}; differ {differ:.3%} of {union}")
    reason = Counter()
    gaps = []
    for k in only_h:
        if k not in nail:
            reason["HMMER only: nail reported nothing"] += 1
        else:
            reason["HMMER only: nail below GA"] += 1
            gaps.append(models[k[1]][0] - max(a[0] for a in nail[k]))
    print(dict(reason), "nail's shortfall, bits, median", round(statistics.median(gaps), 2) if gaps else None)

    same_pair = same_coords = 0
    jaccard, score_diff = [], []
    for k in both:
        pair_same = pair_coords = True
        for r in hmmer[k]:
            cr = columns(r["hmm_from"], r["ali_from"], r["cigar"])
            best, best_j = None, -1.0
            for a in called[k]:
                ca = columns(a[2], a[4], a[6])
                j = len(cr & ca) / max(1, len(cr | ca))
                if j > best_j:
                    best, best_j = a, j
            jaccard.append(best_j)
            score_diff.append(best[0] - r["dom_score"])
            pair_same &= (best[2], best[4], best[6]) == (r["hmm_from"], r["ali_from"], r["cigar"])
            pair_coords &= (best[2], best[3], best[4], best[5]) == (r["hmm_from"], r["hmm_to"], r["ali_from"], r["ali_to"])
        pair_same &= len(called[k]) == len(hmmer[k])
        same_pair += pair_same
        same_coords += pair_coords
    n = max(1, len(both))
    print(f"pairs called by both: identical alignments {same_pair / n:.2%}, same coordinates {same_coords / n:.2%}")
    jaccard.sort()
    score_diff.sort()
    q = lambda v, f: round(v[int(f * (len(v) - 1))], 3) if v else None  # noqa: E731
    print("aligned columns shared per domain (Jaccard): mean", round(sum(jaccard) / max(1, len(jaccard)), 4), "1%", q(jaccard, 0.01), "5%", q(jaccard, 0.05), "median", q(jaccard, 0.5))
    print("nail - HMMER score, bits: 1%", q(score_diff, 0.01), "5%", q(score_diff, 0.05), "median", q(score_diff, 0.5), "95%", q(score_diff, 0.95), "99%", q(score_diff, 0.99))

    # Nail's calls in domtbl columns, so the pipeline's parser resolves overlaps and counts as usual
    nail_dom = out / "nail_ga.domtbl"
    with open(nail_dom, "w") as f:
        for (orf, pf), alis in called.items():
            ga = models[pf]
            for i, a in enumerate(alis, 1):
                f.write(f"{orf} - 0 {pf} {ga[3]} {ga[2]} {a[1]} {a[0]} 0 {i} {len(alis)} {a[1]} {a[1]} {a[0]} 0 {a[2]} {a[3]} {a[4]} {a[5]} {a[4]} {a[5]} 0 -\n")
    a, b = profile(domtbl, parser, out / "hmmer_profile"), profile(nail_dom, parser, out / "nail_profile")
    fams = a.keys() | b.keys()
    moved = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in fams)
    diff = sorted(((b.get(k, 0) - a.get(k, 0), k, a.get(k, 0), b.get(k, 0)) for k in fams if a.get(k, 0) != b.get(k, 0)), key=lambda x: -abs(x[0]))
    print(f"profile: families HMMER {len(a)}, nail {len(b)}, differing {len(diff)}; reads total HMMER {sum(a.values())}, nail {sum(b.values())}; absolute difference {moved} ({moved / max(1, sum(a.values())):.2%})")
    print("largest differences (family, HMMER, nail):", [(k, x, y) for _, k, x, y in diff[:10]])
    with open(out / "profile_differences.tsv", "w") as f:
        f.write("family\thmmer\tnail\n")
        f.writelines(f"{k}\t{x}\t{y}\n" for _, k, x, y in diff)
    verdict = differ < 0.001 and same_pair / n >= 0.99 and not diff
    print("nail alone passes the plan's rule:", verdict)


if __name__ == "__main__":
    main()
