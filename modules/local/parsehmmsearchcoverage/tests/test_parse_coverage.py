"""The parser gives what resolving the whole table at once gives, in the same order.
Run: python modules/local/parsehmmsearchcoverage/tests/test_parse_coverage.py"""

import json
import random
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "bin/hmmer_domtbl_parse_coverage.py"


def expected(hits):
    """Positions on the read and model as sets, every read in memory."""
    by_read = defaultdict(list)
    for hit in hits:
        by_read[hit[0]].append(hit)
    reads, depth = defaultdict(set), {}
    for read, vs in by_read.items():
        taken = set()
        for _, frame, begin, end, model, length, evalue, hmm_from, hmm_to, ali_from, ali_to in sorted(vs, key=lambda h: float(h[6])):
            # residue ali_from to ali_to + 1, as 3 nt each from the frame's end of the read
            at = [(begin + 3 * (x - 1) + frame) if frame > 0 else (end - 3 * (x - 1) - frame) for x in (ali_from, ali_to + 1)]
            nt = set(range(min(at), max(at)))
            if taken & nt:
                continue
            taken |= nt
            reads[model].add(read)
            cover = depth.setdefault(model, [0] * length)
            for i in range(hmm_from - 1, hmm_to):
                cover[i] += 1
    rows = [(m, len(reads[m]), sum(d) / len(d), sum(v > 0 for v in d) / len(d)) for m, d in depth.items()]
    stats = {"reads_mapped": len(by_read), "hmm_count": len(reads), "read_hit_count": sum(map(len, reads.values()))}
    return sorted(rows, key=lambda r: -r[2]), stats


def random_hits(rng, n):
    lengths = [rng.randrange(20, 200) for _ in range(30)]
    for _ in range(n):
        model = rng.randrange(30)
        hmm_from = rng.randrange(1, lengths[model] + 1)
        ali_from = rng.randrange(1, 30)
        yield (
            f"r_{rng.randrange(n // 4)}_length=150", rng.choice([-3, -2, -1, 1, 2, 3]), rng.randrange(1, 4), rng.randrange(120, 151),
            f"PF{model:05}.1", lengths[model], rng.choice(["1e-05", "2.3e-09", "0.001", "0"]),
            hmm_from, rng.randrange(hmm_from, lengths[model] + 1), ali_from, ali_from + rng.randrange(20),
        )  # fmt: skip


def test_parser():
    rng = random.Random(1)
    for n in (4, 400, 4000):
        hits = list(random_hits(rng, n))
        table = "# header\n" + "".join(
            f"{r}_frame={f}_begin={b}_end={e} - 50 name {m} {length} {ev} 30.1 0.0 1 1 1e-9 1e-9 30.0 0.0 {hf} {ht} {af} {at} {af} {at} 0.9 a description\n"
            for r, f, b, e, m, length, ev, hf, ht, af, at in hits
        )
        with tempfile.TemporaryDirectory() as d:
            subprocess.run([sys.executable, SCRIPT, "-i", "-", "-o", "o.tsv", "-s", "o.json"], input=table.encode(), cwd=d, check=True)
            rows, stats = expected(hits)
            want = "# function\tread_count\tcoverage_depth\tcoverage_breadth\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows)
            assert (Path(d) / "o.tsv").read_text() == want
            assert json.loads((Path(d) / "o.json").read_text()) == stats
            assert {p.name for p in Path(d).iterdir()} == {"o.tsv", "o.json"}, "temporary files left"


if __name__ == "__main__":
    test_parser()
    print("ok")
