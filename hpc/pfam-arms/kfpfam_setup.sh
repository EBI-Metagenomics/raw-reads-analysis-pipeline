#!/usr/bin/env bash
#SBATCH --job-name=kfpfam-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=kfpfam-setup.%j.log
#
# One-off setup for the kfpfam arm: the image, the Pfam prep directory and kfpfam's top index, the standard
# mode for reads. Run it on the HPC before run.sh, from this directory, after lookup_setup.sh. Everything
# lands in kfpfam/ beside this script, with the Pfam-A.hmm it was built from; nothing here is part of
# the benchmark's measured cost.
#
#   KFPFAM=~/src/kfpfam KFP=~/src/kmer-functional-profiler sbatch kfpfam_setup.sh
#
#   KFPFAM   a kfpfam checkout, branch prep-any-hmm (its docker/filter.Dockerfile, if there is no image tar)
#   KFP      a kmer-functional-profiler checkout, branch pfam-anchor: top is built by hm-annotate's tables.py
#            and seeds.py, as its PREP does
#   PFAM_HMM the Pfam-A.hmm the kfpfam arm searches: by default the lookup bundle's own copy
#            (lookup/bundle/hmm/Pfam-A.hmm, Pfam 38.1), so the reads' alignments are the bundle's Pfam.
#            A prep and index built from any other file are moved aside and rebuilt.
#   N_MODELS optional: how many models PFAM_HMM must hold
#
# Top: seeds of shape 10110001 at theta 9 (a per-model threshold), a gapped bar per model in
# prep/gapped.txt (KFPFAM_HMMER uses it when it is there), and the pipeline's kfpfam defaults (bar 15,
# gate 14.5, HMMER at 21.25). The same build as hm-annotate's PREP for --kfpfam_mode top.
set -euo pipefail
: "${KFPFAM:?set KFPFAM to a kfpfam checkout}"
: "${KFP:?set KFP to a kmer-functional-profiler checkout}"
PFAM_HMM=${PFAM_HMM:-$PWD/lookup/bundle/hmm/Pfam-A.hmm}
[[ -s $PFAM_HMM ]] || { echo "no $PFAM_HMM: run lookup_setup.sh first, or set PFAM_HMM" >&2; exit 1; }
PFAM_HMM=$(realpath "$PFAM_HMM")
HA=$(realpath "$KFP")/workflows/hmm-annotate
[[ -s $HA/seeds.py ]] || { echo "no $HA/seeds.py: KFP is not a kmer-functional-profiler checkout" >&2; exit 1; }
OUT=$PWD/kfpfam
SIF=$OUT/kfpfam-filter.sif
mkdir -p "$OUT"

# kfpfam is published nowhere, so there is nothing to pull: the image is built where Docker is,
# then copied here. Codon's nodes have no Docker and `singularity build` needs root or fakeroot.
# Rebuilt when the tar is newer, so a re-exported image doesn't leave a stale .sif behind.
if [[ ! -s $SIF || ( -s $OUT/kfpfam-filter.tar && $OUT/kfpfam-filter.tar -nt $SIF ) ]]; then
    rm -f "$SIF"
    version=$(git -C "$KFPFAM" rev-parse --short HEAD 2>/dev/null || echo unknown)
    if [[ -s $OUT/kfpfam-filter.tar ]]; then
        # Exported with `docker save` elsewhere: convert it here, no Docker needed.
        singularity build "$SIF" "docker-archive://$OUT/kfpfam-filter.tar"
    else
        cat >&2 <<MSG
no image. Build it where Docker is, export it, and copy the tar to $OUT:
  docker build -f docker/filter.Dockerfile --build-arg KFPFAM_VERSION=$version \\
      --platform linux/amd64 -t kfpfam-filter $KFPFAM
  docker save kfpfam-filter:latest -o kfpfam-filter.tar
Then re-run this script, which converts it to $SIF.
The image must be amd64: Codon is x86_64, and the filter picks its SIMD path at run time.
MSG
        exit 1
    fi
fi

run () { singularity exec --bind "$PWD:$PWD" --bind "$HA:$HA" --pwd "$OUT" "$SIF" "$@"; }

# The pipeline joins seeds to models by name, so prep must hold PFAM_HMM's models, in its order. The
# first Codon build was from 39.0 while the pipeline searched 38.2, and 183 families were lost silently.
zcat -f "$PFAM_HMM" | awk '$1 == "NAME" { n = $2 } $1 == "ACC" { split($2, a, "."); print a[1] "\t" n }' > "$OUT/models.expected"
n=$(wc -l < "$OUT/models.expected")
if [[ -n ${N_MODELS:-} ]] && (( n != N_MODELS )); then
    echo "$PFAM_HMM holds $n models, not $N_MODELS: it is another Pfam release" >&2
    exit 1
fi
# names can agree between releases whose profiles don't, so the file itself is compared too
sum=$(zcat -f "$PFAM_HMM" | sha256sum | cut -d' ' -f1)
if [[ -s $OUT/prep/models.tsv && $(cat "$OUT/prep/hmm.sha256" 2>/dev/null) != "$sum" ]]; then
    old=$OUT/stale.$(date +%Y%m%d-%H%M%S)
    echo "prep and index were not built from $PFAM_HMM: moving them to $old and rebuilding"
    mkdir "$old"
    mv "$OUT/prep" "$old/"
    [[ ! -d $OUT/index ]] || mv "$OUT/index" "$old/"
fi

if [[ ! -s $OUT/index/top.tidx ]]; then
    rm -rf "$OUT/prep" "$OUT/index" && mkdir -p "$OUT/index"
    export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
    # kept in kfpfam/, which the containers see: run.sh's kfpfam arm searches this copy, the index's own
    zcat -f "$PFAM_HMM" > "$OUT/Pfam-A.hmm"
    # what hm-annotate's PREP does for --kfpfam_mode top
    run python3 /opt/kfpfam/prep.py Pfam-A.hmm - prep
    if ! cmp -s "$OUT/models.expected" <(tail -n +2 "$OUT/prep/models.tsv" | cut -f1,2); then
        echo "prep.py's models are not $PFAM_HMM's" >&2
        exit 1
    fi
    run python3 "$HA/tables.py" thresholds --hmm Pfam-A.hmm --models prep/models.tsv --cutoff gathering --out thresholds.tsv
    run python3 "$HA/seeds.py" bars --models prep/models.tsv --thresholds thresholds.tsv --ga ga.tsv --gapped prep/gapped.txt --slack 0
    run filter seeds prep ga.tsv seeds.f32 10110001
    run python3 "$HA/seeds.py" choose --table seeds.f32 --shape 10110001 --theta 9 --out T.txt
    run env FILTER_T_FILE=T.txt filter build prep 10110001 0 index/top.idx
    run filter tile prep index/top.idx index/top.tidx.part
    mv "$OUT/index/top.tidx.part" "$OUT/index/top.tidx"
    rm -f "$OUT/index/top.idx" "$OUT/seeds.f32"
    echo "$sum" > "$OUT/prep/hmm.sha256"   # last, so a build that failed is not taken for a good one
fi

echo "Pfam:      $OUT/Pfam-A.hmm (from $PFAM_HMM), $n models"
echo "prep:      $OUT/prep"
echo "index:     $OUT/index/top.tidx"
echo "image:     $SIF"
