#!/usr/bin/env bash
# PFAM_LOOKUP in two parts and PFAM_LOOKUP_MERGE, not their stubs, on a fixture bundle (kmer-functional-profiler's lookup-index test
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
include { PFAM_LOOKUP_MERGE } from '$root/modules/local/pfam_lookup_merge/main'
workflow {
    def fx = '$root/tests/fixtures/pfam_lookup'
    def parts = channel.of(1, 2).map { i -> tuple([id: 'test', lookup_part: i], file("\$fx/alignments.parquet"), "\${i}/2") }
    PFAM_LOOKUP(parts, file("\$fx/bundle"), 0.4, 2000, 1000, 40)
    PFAM_LOOKUP_MERGE(PFAM_LOOKUP.out.reads.map { _m, r -> r }.collect().map { r -> tuple([id: 'test'], r) }, file("\$fx/bundle"), 0.4)
}
NF
cat > "$tmp/nextflow.config" <<CFG
process.publishDir = '$tmp/out'
process { withName: PFAM_LOOKUP { ext.prefix = { "\${meta.id}_pfam_lookup.part\${meta.lookup_part}" } } }
process { withName: PFAM_LOOKUP_MERGE { publishDir = [path: '$tmp/out', saveAs: { f -> f == 'versions.yml' ? null : f }] } }
process.cpus = 1
process.memory = 1.GB
process.ext.bundle_py = '$BUNDLE_PY'
env.LOOKUP_BIN = '$LOOKUP_BIN'
CFG
(cd "$tmp" && PATH="$PYTHON_BIN:$PATH" nextflow run main.nf -ansi-log false > nf.log) || { cat "$tmp/nf.log"; exit 1; }
grep -q '"reads_called": 18' "$tmp/out/test_pfam_lookup.summary.json"
grep -q '^K00001+K00002	6	' "$tmp/out/test_pfam_lookup.ko.tsv"
test -s "$tmp/out/test_pfam_lookup.reads.parquet"
test -s "$tmp/out/test_pfam_lookup.part1.reads.parquet" && test -s "$tmp/out/test_pfam_lookup.part2.reads.parquet"
grep -q 'format 5, Pfam 38.2' "$tmp/out/versions.yml"
echo "PFAM_LOOKUP's outputs are as expected"
