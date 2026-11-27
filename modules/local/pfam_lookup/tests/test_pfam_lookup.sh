#!/usr/bin/env bash
# PFAM_LOOKUP itself, not its stub, on a fixture bundle (kmer-functional-profiler's lookup-index test
# output), outside its image. Run from the repository root:
#   BUNDLE_PY=<kmer-functional-profiler>/workflows/pfam-anchor/lookup/bundle.py LOOKUP_BIN=<lookup> \
#       PYTHON_BIN=<directory of a python3 with polars, pyarrow and numpy> bash modules/local/pfam_lookup/tests/test_pfam_lookup.sh
set -euo pipefail
: "${BUNDLE_PY:?}" "${LOOKUP_BIN:?}" "${PYTHON_BIN:?}"
root=$PWD
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cat > "$tmp/main.nf" <<NF
params.pfam_lookup_sif = ''
include { PFAM_LOOKUP } from '$root/modules/local/pfam_lookup/main'
workflow {
    def fx = '$root/tests/fixtures/pfam_lookup'
    PFAM_LOOKUP([[id: 'test'], file("\$fx/alignments.parquet")], file("\$fx/bundle"), 0.4, 2000, 1000, 40)
}
NF
cat > "$tmp/nextflow.config" <<CFG
process.publishDir = '$tmp/out'
process.cpus = 1
process.memory = 1.GB
process.ext.bundle_py = '$BUNDLE_PY'
env.LOOKUP_BIN = '$LOOKUP_BIN'
CFG
(cd "$tmp" && PATH="$PYTHON_BIN:$PATH" nextflow run main.nf -ansi-log false > nf.log) || { cat "$tmp/nf.log"; exit 1; }
grep -q '"reads_called": 18' "$tmp/out/test_pfam_lookup.summary.json"
grep -q '^K00001+K00002	6	' "$tmp/out/test_pfam_lookup.ko.tsv"
test -s "$tmp/out/test_pfam_lookup.reads.parquet"
grep -q 'format 5, Pfam 38.2' "$tmp/out/versions.yml"
echo "PFAM_LOOKUP's outputs are as expected"
