"""merge_parquet.py gives what one sort of everything gives. Needs pyarrow.
Run: python modules/local/mergeparquet/tests/test_merge_parquet.py"""

import random
import sys
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from merge_parquet import merge, range_starts  # noqa: E402

KEYS = [("pfam", "ascending"), ("read", "ascending")]


def test_merge():
    rng = random.Random(1)
    with tempfile.TemporaryDirectory() as d:
        files, tables = [], []
        for i in range(3):
            # PF00000 alone is larger than a pass
            t = pa.table({"read": [f"r{i}.{n}" for n in range(200)],
                          "pfam": [f"PF{max(0, rng.randrange(-20, 20)):05}" for _ in range(200)]})
            tables.append(t)
            files.append(f"{d}/{i}.parquet")
            pq.write_table(t, files[-1], row_group_size=50)
        merge(files, f"{d}/out.parquet", rows=60)
        assert len(range_starts(files, 60)) > 3
        assert not [x for x in Path(d).iterdir() if x.is_dir()], "temporary files left"
        assert pq.read_table(f"{d}/out.parquet").equals(pa.concat_tables(tables).combine_chunks().sort_by(KEYS))
        pq.write_table(tables[0].slice(0, 0), files[0])
        merge(files[:1], f"{d}/empty.parquet")
        assert pq.read_table(f"{d}/empty.parquet").equals(tables[0].slice(0, 0))


if __name__ == "__main__":
    test_merge()
    print("ok")
