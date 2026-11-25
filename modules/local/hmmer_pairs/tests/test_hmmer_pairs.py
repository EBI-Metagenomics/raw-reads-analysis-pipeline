"""hmmer_pairs.py's Parquet output. Needs pyhmmer and pyarrow.
Run: python modules/local/hmmer_pairs/tests/test_hmmer_pairs.py"""

import gzip
import io
import re
import struct
import sys
import tempfile
from array import array
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "modules/local/hmmer_pairs/bin"))
import hmmer_pairs as hp  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402


FX = ROOT / "tests/fixtures"
ORFS, HMM = FX / "hmmer_pairs/orfs.faa.gz", str(FX / "dbs/pfam_mini/Pfam-A.hmm")


def run(jobs, out, threads, parquet=True, weak=0):
    path = out / "x.parquet"
    with hp.Output(out / "x.domtbl.gz", path if parquet else None) as output, ThreadPoolExecutor(threads) as pool:
        for j in jobs:
            hp.run(j, 34, threads, pool, output, weak)
    return pq.read_table(path) if parquet else None


def run_on_fixture(out, threads, parquet=True, weak=0):
    pairs, ids = hp.read_seeds(FX / "hmmer_pairs/seeds.tsv")
    seqs = hp.read_orfs(ORFS, ids)
    return run([hp.seeded(pairs, seqs, HMM)], out, threads, parquet, weak)


def stream(block):
    """The fixture's seeds as `filter join` writes them, `block` ORFs at a time, and models.tsv's text."""
    with gzip.open(ORFS, "rt") as f:
        orfs = {r.split()[0]: len("".join(r.split("\n")[1:])) for r in f.read()[1:].split("\n>")}
    with open(HMM) as f:
        models = [line.split()[1] for line in f if line.startswith("NAME ")]
    with open(FX / "hmmer_pairs/seeds.tsv") as f:
        seeds = {tuple(line.split()[:2]) for line in f}
    names, data = list(orfs), b""
    for start in range(0, len(names), block):
        at = {n: i for i, n in enumerate(names[start : start + block])}
        pairs = sorted((models.index(m), at[o]) for m, o in seeds if o in at)
        flat = array("I", [x for m, o in pairs for x in (o, m)]).tobytes()
        data += struct.pack("<3Q", len(at), sum(orfs[n] for n in at), len(flat)) + flat
    return data, "acc\tname\n" + "".join(f"PF\t{m}\n" for m in models)


def run_stream(data, models, out, size):
    (out / "models.tsv").write_text(models)
    optimized = hp.profiles([HMM], out / "models.tsv")
    with hp.open_orfs(ORFS) as orfs:
        blocks = hp.read_blocks(io.BytesIO(data), orfs)
        return run(([(optimized[m], t[m]) for m in sorted(t)] for t in hp.batches(blocks, size)), out, 2)


def domtbl(out):
    """Rows but for the c-Evalue, which depends on how many ORFs were searched together."""
    with gzip.open(out / "x.domtbl.gz", "rt") as f:
        return sorted(tuple(c for i, c in enumerate(line.split()) if i != 11) for line in f)


def test_stream_gives_what_seeds_give():
    keys = [(k, "ascending") for k in ("pfam", "read", "orf_begin", "frame", "domain")]
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        a, b = Path(a), Path(b)
        want = run_on_fixture(a, 2).sort_by(keys)
        assert want.num_rows > 0
        for block, size in ((100, 1), (100, 10**6), (10**6, 1)):
            assert run_stream(*stream(block), b, size).sort_by(keys).equals(want), (block, size)
            assert domtbl(a) == domtbl(b), (block, size)


def test_a_bad_stream_is_an_error():
    """A filter that dies, or that read other ORFs or another Pfam, must not give a smaller profile."""
    data, models = stream(100)
    first = 24 + struct.unpack("<3Q", data[:24])[2]
    wrong = struct.pack("<3Q", 100, 1, 0)
    cases = ((data[:-1], models), (data[:first], models), (data[: first + 5], models), (b"", models),
             (wrong + data[first:], models), (data, models.replace("AAA\n", "BBB\n", 1)))
    for i, (d, m) in enumerate(cases):
        with tempfile.TemporaryDirectory() as out:
            try:
                run_stream(d, m, Path(out), 10**6)
            except ValueError:
                continue
            raise AssertionError(f"case {i} was not refused")


def test_alignments_are_optional():
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        run_on_fixture(Path(a), 2)
        run_on_fixture(Path(b), 2, parquet=False)
        assert domtbl(Path(a)) == domtbl(Path(b))
        assert not (Path(b) / "x.parquet").exists()


def test_weak_alignments():
    """They are rows added under the gathering thresholds: the domtbl is the same byte for byte, and so are the other rows."""
    key = ["pfam", "read", "frame", "orf_begin", "env_from", "env_to"]
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b, tempfile.TemporaryDirectory() as c:
        ga = run_on_fixture(Path(a), 2).to_pylist()
        loose = run_on_fixture(Path(b), 2, weak=1000.0).to_pylist()
        tight = run_on_fixture(Path(c), 1, weak=1e-3).to_pylist()
        for x in (b, c):
            assert gzip.open(Path(a) / "x.domtbl.gz").read() == gzip.open(Path(x) / "x.domtbl.gz").read()
    assert ga and all(r["above_ga"] for r in ga)
    for rows, limit in ((loose, 1000.0), (tight, 1e-3)):
        assert [r for r in rows if r["above_ga"]] == ga
        weak = [r for r in rows if not r["above_ga"]]
        assert all(r["seq_evalue"] <= limit and r["i_evalue"] <= limit for r in weak)
        assert len({tuple(r[k] for k in key) for r in rows}) == len(rows), "a domain is there twice"
        assert len({(r["pfam"], r["read"], r["frame"], r["orf_begin"], r["domain"]) for r in rows}) == len(rows)
    assert any(not r["above_ga"] for r in loose)


def test_the_models_version():
    """Columns only mean something for one version of a model, so each row says which."""
    with open(HMM) as f:
        versions = dict(line.split()[1].split(".") for line in f if line.startswith("ACC "))
    with tempfile.TemporaryDirectory() as d:
        for r in run_on_fixture(Path(d), 2).to_pylist():
            assert r["pfam_version"] == int(versions[r["pfam"]])


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


def test_seeds_from_another_release_are_an_error():
    """A name the whole database lacks is an error; one a chunk lacks is another chunk's."""
    with tempfile.TemporaryDirectory() as d:
        for whole in (True, False):
            try:
                run([hp.seeded({"Not_in_this_Pfam": [0]}, [None], HMM, whole)], Path(d), 1)
                assert not whole, "an unknown model was skipped"
            except ValueError as e:
                assert whole and "Not_in_this_Pfam" in str(e), e


if __name__ == "__main__":
    test_schema_and_alignment_lengths()
    test_read_span_covers_the_right_codons()
    test_output_does_not_depend_on_threads()
    test_seeds_from_another_release_are_an_error()
    test_stream_gives_what_seeds_give()
    test_a_bad_stream_is_an_error()
    test_alignments_are_optional()
    test_weak_alignments()
    test_the_models_version()
    print("ok")
