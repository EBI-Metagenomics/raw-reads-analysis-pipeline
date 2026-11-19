"""select_frames.py's frame logic. Run: python3 modules/local/selectframes/tests/test_select_frames.py"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "modules/local/selectframes/bin"))
import select_frames as sf  # noqa: E402

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def test_nucleotides_translate_to_seqkit_orfs():
    """Both strands, on real seqkit output: the coordinates give back the ORF's own protein."""
    fx = ROOT / "tests/fixtures"
    n = {False: 0, True: 0}
    for read, seq, rows in sf.by_read(fx / "profile_hmmsearch_pfam/ERR10889147.mini.fasta.gz", fx / "hmmer_pairs/orfs.faa.gz"):
        for header, protein in rows:
            name, frame, nt = sf.nucleotides(seq, header)
            assert name == read
            got = "".join(CODON.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3))
            assert got == protein, header
            n[frame[0]] += 1
    assert n[False] > 100 and n[True] > 100, n


def test_call_frames_are_orf_frames():
    """A call in frame with an ORF has the same key on both strands."""
    path = Path(sys.argv[0]).with_name("calls_test.out")
    path.write_text(">r\n4\t60\t+\t1\t1.0\tI:\tD:\n>s\n3\t62\t-\t1\t1.0\tI:\tD:\n>t\n1\t30\t+\t1\t1.0\tI:\tD:\n32\t60\t+\t2\t1.0\tI:\tD:\n")
    try:
        calls = sf.read_calls(path)
    finally:
        path.unlink()
    assert calls["r"] == {(False, 0)}
    seq = "A" * 70
    assert sf.nucleotides(seq, "r_frame=1_begin=4_end=60")[1] in calls["r"]
    assert calls["s"] == {(True, 2)}
    assert sf.nucleotides(seq, "s_frame=-1_begin=3_end=62")[1] in calls["s"]
    assert calls["t"] == {(False, 0), (False, 1)}


if __name__ == "__main__":
    test_nucleotides_translate_to_seqkit_orfs()
    test_call_frames_are_orf_frames()
    print("ok")
