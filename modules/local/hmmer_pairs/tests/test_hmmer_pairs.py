"""hmmer_pairs.py's Parquet output. Needs pyhmmer and pyarrow.
Run: python modules/local/hmmer_pairs/tests/test_hmmer_pairs.py"""

import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "modules/local/hmmer_pairs/bin"))
import hmmer_pairs as hp  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402


def run_on_fixture(out, threads):
    fx = ROOT / "tests/fixtures"
    pairs, ids = hp.read_seeds(fx / "hmmer_pairs/seeds.tsv")
    seqs = hp.read_orfs(fx / "hmmer_pairs/orfs.faa.gz", ids)
    hp.run(pairs, seqs, str(fx / "dbs/pfam_mini/Pfam-A.hmm"), 34, threads, str(out / "x.domtbl.gz"), str(out / "x.parquet"))
    return pq.read_table(out / "x.parquet")


def test_schema_and_alignment_lengths():
    with tempfile.TemporaryDirectory() as d:
        table = run_on_fixture(Path(d), 2)
    assert table.schema.equals(hp.SCHEMA), table.schema
    assert table.num_rows > 0
    for r in table.to_pylist():
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MID])", r["cigar"])]
        assert sum(n for n, op in ops if op != "D") == len(r["residues"]) == r["ali_to"] - r["ali_from"] + 1
        assert sum(n for n, op in ops if op != "I") == r["hmm_to"] - r["hmm_from"] + 1
        assert sum(n for n, _ in ops) == len(r["pp"])
    assert table.column("pfam").to_pylist() == sorted(table.column("pfam").to_pylist())


def test_output_does_not_depend_on_threads():
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        assert run_on_fixture(Path(a), 1).equals(run_on_fixture(Path(b), 4))


if __name__ == "__main__":
    test_schema_and_alignment_lengths()
    test_output_does_not_depend_on_threads()
    print("ok")
