#!/usr/bin/env bash
#SBATCH --job-name=weefams-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=weefams-setup.%j.log
#
# One-off setup for the weefam arms: kfpfam's top mode on mgnifams-explore's weefams. Run it on the HPC after
# kfpfam_setup.sh (whose image it uses), from this directory. Everything lands in weefams/ beside this script.
#
#   WEEFAMS=<runs>/weefams/full/results/wcollect/weefams.hmm KFP=~/src/kmer-functional-profiler ./weefams_setup.sh
#
#   WEEFAMS  the weefam library: the accepted weefams, each with its cut as GA
#   KFP      a kmer-functional-profiler checkout: top is built by hm-annotate's tables.py and seeds.py, as its PREP does
#
# Writes weefams/weefams.hmm (the library with ACC lines: the pipeline counts reads by accession, and a weefam has
# none, so its NAME is its ACC), weefams/prep with gapped.txt (a gapped bar per model, as top's) and
# weefams/index/top.tidx (top's seeds: theta 9, spaced seed 10110001). weefams.config points the pipeline at them.
set -euo pipefail
: "${WEEFAMS:?set WEEFAMS to weefams.hmm}"
: "${KFP:?set KFP to a kmer-functional-profiler checkout}"
OUT=$PWD/weefams
SIF=$PWD/kfpfam/kfpfam-filter.sif
HA=$KFP/workflows/hmm-annotate
[[ -s $SIF ]] || { echo "no $SIF: run kfpfam_setup.sh first" >&2; exit 1; }
mkdir -p "$OUT"
run () { singularity exec --bind "$PWD:$PWD" --bind "$HA:$HA" --pwd "$OUT" "$SIF" "$@"; }

# A new library replaces the old prep and index
sum=$(zcat -f "$WEEFAMS" | sha256sum | cut -d' ' -f1)
if [[ -s $OUT/prep/models.tsv && $(cat "$OUT/source.sha256" 2>/dev/null) != "$sum" ]]; then
    old=$OUT/stale.$(date +%Y%m%d-%H%M%S)
    echo "weefams/ was not built from $WEEFAMS: moving it to $old"
    mkdir "$old"
    mv "$OUT/prep" "$OUT/index" "$OUT/weefams.hmm" "$old/" 2>/dev/null || true
fi

if [[ ! -s $OUT/index/top.tidx ]]; then
    zcat -f "$WEEFAMS" | awk '$1 == "ACC" { next } { print } $1 == "NAME" { print "ACC   " $2 }' > "$OUT/weefams.hmm.part"
    if [[ $(grep -c '^NAME ' "$OUT/weefams.hmm.part") != $(grep '^NAME ' "$OUT/weefams.hmm.part" | sort -u | wc -l) ]]; then
        echo "$WEEFAMS: two weefams share a name" >&2
        exit 1
    fi
    mv "$OUT/weefams.hmm.part" "$OUT/weefams.hmm"
    rm -rf "$OUT/prep" "$OUT/index" && mkdir "$OUT/index"
    export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
    # what hm-annotate's PREP does for --kfpfam_mode top
    run python3 /opt/kfpfam/prep.py weefams.hmm - prep
    run python3 "$HA/tables.py" thresholds --hmm weefams.hmm --models prep/models.tsv --cutoff gathering --out thresholds.tsv
    run python3 "$HA/seeds.py" bars --models prep/models.tsv --thresholds thresholds.tsv --ga ga.tsv --gapped prep/gapped.txt --slack 0
    run filter seeds prep ga.tsv seeds.f32 10110001
    run python3 "$HA/seeds.py" choose --table seeds.f32 --shape 10110001 --theta 9 --out T.txt
    run env FILTER_T_FILE=T.txt filter build prep 10110001 0 index/top.idx
    run filter tile prep index/top.idx index/top.tidx.part
    mv "$OUT/index/top.tidx.part" "$OUT/index/top.tidx"
    rm -f "$OUT/index/top.idx" "$OUT/seeds.f32"
    echo "$sum" > "$OUT/source.sha256"   # last, so a build that failed is not taken for a good one
fi

echo "weefams:   $OUT/weefams.hmm, $(grep -c '^NAME ' "$OUT/weefams.hmm") models"
echo "prep:      $OUT/prep"
echo "index:     $OUT/index/top.tidx"
