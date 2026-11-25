#!/usr/bin/env bash
#SBATCH --job-name=pfam-arms
#SBATCH --output=%x.%j.log
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# V5 of docs/plans/pfam-orf-filter-nail.md: four arms of Pfam profiling on samplesheet.csv.
# Submit from this directory: cd hpc/pfam-arms && sbatch run.sh
# Re-submitting resumes. Each arm reuses the others' QC, decontamination and translation.
# Read:
#   results/compare.tsv   per sample and arm, against `current`: read and family differences, worst families
#   results/families.tsv  reads per family in each arm
#   results/cost.tsv      core-hours and peak memory per Pfam step
#   results/<arm>/<sample>/function-summary/pfam/  each arm's published files, with filter stats
set -euo pipefail
[[ -f run.sh && -f site.config && -f unchunked_pfam.config ]] || { echo "submit from hpc/pfam-arms" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g' NXF_VER=25.10.4
PIPELINE=../..

# The nail arms search one whole Pfam-A.hmm, where the hmmsearch arms use the site's chunks.
# Set PFAM_HMM to a single Pfam-A.hmm if you have one; otherwise set PFAM_CHUNKS to the chunk
# files the site config serves (a glob is fine) and this builds it once, beside this script.
# An HMM file is its models one after another, so concatenating the chunks gives the database
# the hmmsearch arms search, and -Z below is counted from the result rather than assumed.
PFAM_HMM=${PFAM_HMM:-$PWD/pfam_unchunked/Pfam-A.hmm}
if [[ ! -s $PFAM_HMM ]]; then
    [[ -n ${PFAM_CHUNKS:-} ]] || {
        echo "set PFAM_HMM to one Pfam-A.hmm, or PFAM_CHUNKS to the chunk HMMs to join" >&2
        exit 1
    }
    mkdir -p "$(dirname "$PFAM_HMM")"
    cat $PFAM_CHUNKS > "$PFAM_HMM.part"
    mv "$PFAM_HMM.part" "$PFAM_HMM"
fi
export PFAM_HMM_DIR=$(dirname "$PFAM_HMM")
export PFAM_HMM_NAME=$(basename "$PFAM_HMM")
export PFAM_NUM_MODELS=$(grep -c '^NAME ' "$PFAM_HMM")
echo "nail arms search $PFAM_HMM, $PFAM_NUM_MODELS models"

arm() {
    local name=$1
    shift
    nextflow run "$PIPELINE" -profile singularity -c site.config -resume \
        --samplesheet samplesheet.csv --outdir "results/$name" \
        --skip_functional false --skip_taxonomic true \
        --hmmsearch_subsampling=-1 \
        --trace_report_suffix arm --singularity_cachedir singularity_cache "$@"
}

# Seeded arms first: they are the cheap ones, and whichever runs first pays for the shared QC,
# decontamination and merge, so the hmmsearch arms start from a warm cache rather than the other
# way round. With subsampling off, decide from the seeded arms' real cost whether the hmmsearch
# arms are worth running at full depth at all. kfpfam goes first of all: its prefilter is the
# cheapest, so a failure there surfaces before anything expensive is queued.

# kfpfam's two modes, on one whole Pfam like the nail arms: the index covers the whole database,
# and chunks would only split HMMER_PAIRS. The index carries its own k-mer bit threshold, so each
# mode has its own and the numbers must match it: A is nail's cost at higher recall, B is nail's
# recall at less cost.
if [[ -s kfpfam/index/k4T6.5.tidx ]]; then
    kf="--kfpfam_sif $PWD/kfpfam/kfpfam-filter.sif --kfpfam_prep $PWD/kfpfam/prep"
    # The MGnify lookup on these arms' alignments, when its image (lookup_setup.sh) and bundle are there. The
    # bundle must be of the Pfam searched here, or PFAM_LOOKUP stops. Everything upstream is cached.
    if [[ -s lookup/bundle/manifest.json && -s lookup/kfp-lookup.sif ]]; then
        kf+=" --pfam_lookup true --pfam_lookup_bundle $(realpath lookup/bundle) --pfam_lookup_sif $PWD/lookup/kfp-lookup.sif"
    else
        echo "no lookup/bundle and lookup/kfp-lookup.sif: the kfpfam arms run without the MGnify lookup" >&2
    fi
    arm kfpfam_a -c unchunked_pfam.config $kf --pfam_search_engine kfpfam \
        --kfpfam_index "$PWD/kfpfam/index/k4T6.5.tidx" \
        --kfpfam_bar 17 --kfpfam_gate 16.5 --kfpfam_min_gapped 21.25
    arm kfpfam_b -c unchunked_pfam.config $kf --pfam_search_engine kfpfam \
        --kfpfam_index "$PWD/kfpfam/index/k4T7.5.tidx" \
        --kfpfam_bar 20 --kfpfam_gate 16.5 --kfpfam_min_gapped 26.25
else
    echo "no kfpfam index: run kfpfam_setup.sh to add the kfpfam arms" >&2
fi

arm nail -c unchunked_pfam.config --pfam_search_engine nail
# nail at its loose seeding filter, so nail and kfpfam A are compared at matched recall rather
# than at matched defaults. V2 measured about 334 candidate pairs per read here against 65.
arm nail_loose -c unchunked_pfam.config --pfam_search_engine nail --pfam_nail_seed_pvalue 0.1
arm both -c unchunked_pfam.config --pfam_search_engine nail --pfam_orf_filter true

# The expensive pair, last: about 20x B11's measured 550 and 358 core-hours now that every merged
# read is searched rather than a 2 M subsample.
arm current
arm filter --pfam_orf_filter true

# `current` first, so compare_profiles.py measures every other arm against it
ARMS="current filter nail both nail_loose"
if [[ -s kfpfam/index/k4T6.5.tidx ]]; then
    ARMS="$ARMS kfpfam_a kfpfam_b"
fi

: > results/compare.tsv
: > results/families.tsv
for sample in $(tail -n +2 samplesheet.csv | cut -d, -f1); do
    args=()
    for a in $ARMS; do
        args+=("$a=results/$a/$sample/function-summary/pfam/${sample}_pfam.txt.gz")
    done
    python3 "$PIPELINE/bin/compare_profiles.py" "$sample" "${args[@]}" --families "results/$sample.families.tsv" >> results/compare.tsv
    cat "results/$sample.families.tsv" >> results/families.tsv
    rm "results/$sample.families.tsv"
done
cost=()
for a in $ARMS; do cost+=("$a=results/$a/pipeline_info/execution_trace_arm.txt"); done
python3 cost.py "${cost[@]}" > results/cost.tsv
echo finished
