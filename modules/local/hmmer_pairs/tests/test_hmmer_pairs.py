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
        # The read span is the aligned residues' nucleotides, so three per residue, inside the ORF
        assert r["read_to"] - r["read_from"] + 1 == 3 * (r["ali_to"] - r["ali_from"] + 1)
        assert min(r["orf_begin"], r["orf_end"]) <= r["read_from"] <= r["read_to"] <= max(r["orf_begin"], r["orf_end"])
    assert table.column("pfam").to_pylist() == sorted(table.column("pfam").to_pylist())


def test_read_span_covers_the_right_codons():
    """Every residue is three nucleotides, counted up from `begin` or down from `end`."""
    for frame, begin, end in ((1, 1, 123), (2, 2, 124), (3, 3, 125), (-1, 3, 125), (-2, 2, 124), (-3, 1, 123)):
        residues = (end - begin + 1) // 3
        for ali_from, ali_to in ((1, 1), (1, 5), (4, 20), (residues, residues), (1, residues)):
            got = hp.read_span(frame, begin, end, ali_from, ali_to)
            assert got[1] - got[0] + 1 == 3 * (ali_to - ali_from + 1), (frame, ali_from, ali_to, got)
            assert begin <= got[0] <= got[1] <= end, (frame, ali_from, ali_to, got)
        # the whole ORF, and each end of it
        assert hp.read_span(frame, begin, end, 1, residues) == (begin, end)
        first, last = hp.read_span(frame, begin, end, 1, 1), hp.read_span(frame, begin, end, residues, residues)
        assert (first, last) == (((begin, begin + 2), (end - 2, end)) if frame > 0
                                 else ((end - 2, end), (begin, begin + 2)))


def test_output_does_not_depend_on_threads():
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        assert run_on_fixture(Path(a), 1).equals(run_on_fixture(Path(b), 4))


if __name__ == "__main__":
    test_schema_and_alignment_lengths()
    test_read_span_covers_the_right_codons()
    test_output_does_not_depend_on_threads()
    print("ok")
