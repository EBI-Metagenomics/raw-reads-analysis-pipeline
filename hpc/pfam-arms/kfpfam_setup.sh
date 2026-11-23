#!/usr/bin/env bash
#SBATCH --job-name=kfpfam-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=kfpfam-setup.%j.log
#
# One-off setup for the kfpfam arms: the image, the Pfam prep directory and one tiled index per
# mode. Run it on the HPC before run.sh, from this directory. Everything lands in kfpfam/ beside
# this script, and nothing here is part of the benchmark's measured cost.
#
#   KFPFAM=~/src/kfpfam PFAM_HMM=pfam_unchunked/Pfam-A.hmm ./kfpfam_setup.sh
#
#   KFPFAM   a kfpfam checkout (its docker/filter.Dockerfile and measure/prep.py are used)
#   PFAM_HMM the same Pfam-A.hmm the arms search, so model names and indices line up with the
#            pipeline's -Z and with HMMER_PAIRS; run.sh builds one from the site's chunks. A prep
#            and index built from any other file are moved aside and rebuilt.
#   N_MODELS optional: how many models it must hold, the pipeline's num_models (53275, Pfam 39.0)
#   CLANS    optional: Pfam-A.clans.tsv.gz. prep.py uses it for one column of models.tsv, and
#            neither the filter nor this pipeline reads that column, so it defaults to empty.
#
# The index encodes its own k-mer bit threshold, so each mode needs its own: Mode A (nail's cost,
# higher recall) uses 6.5 with bar 17 and HMMER at 21.25, Mode B (nail's recall, less cost) uses
# 7.5 with bar 20 and HMMER at 26.25. run.sh passes those numbers; this builds both indexes.
set -euo pipefail
: "${KFPFAM:?set KFPFAM to a kfpfam checkout}"
: "${PFAM_HMM:?set PFAM_HMM to the Pfam-A.hmm the arms search}"
CLANS=${CLANS:-}
N_MODELS=${N_MODELS:-53275}
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

run () { singularity exec --bind "$PWD:$PWD" --pwd "$PWD" "$SIF" "$@"; }

# The pipeline joins seeds to models by name, so prep must hold PFAM_HMM's models, in its order. The
# first Codon build was from 39.0 while the pipeline searched 38.2, and 183 families were lost silently.
# Build from the file the arms search, not from current_release, which moves.
zcat -f "$PFAM_HMM" | awk '$1 == "NAME" { n = $2 } $1 == "ACC" { split($2, a, "."); print a[1] "\t" n }' > "$OUT/models.expected"
n=$(wc -l < "$OUT/models.expected")
if (( n != N_MODELS )); then
    echo "$PFAM_HMM holds $n models, not the pipeline's $N_MODELS: it is another Pfam release (set N_MODELS if it is the one the arms search)" >&2
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

if [[ ! -s $OUT/prep/models.tsv ]]; then
    hmm=$(realpath "$PFAM_HMM")
    gz=$hmm
    [[ $hmm == *.gz ]] || { gzip -c "$hmm" > "$OUT/Pfam-A.hmm.gz"; gz=$OUT/Pfam-A.hmm.gz; }
    if [[ -z $CLANS ]]; then
        # prep.py wants the file, but only to fill models.tsv's `clan` column: the filter loads
        # profile.f32 and offs.u32 alone, and the pipeline reads the `name` column.
        : | gzip -c > "$OUT/clans.tsv.gz"
        CLANS=$OUT/clans.tsv.gz
    fi
    run python3 /opt/kfpfam/prep.py "$gz" "$(realpath "$CLANS")" "$OUT/prep"
    if ! cmp -s "$OUT/models.expected" <(tail -n +2 "$OUT/prep/models.tsv" | cut -f1,2); then
        echo "prep.py's models are not $PFAM_HMM's" >&2
        exit 1
    fi
    echo "$sum" > "$OUT/prep/hmm.sha256"   # last, so a build that failed is not taken for a good one
fi

for spec in "6.5 k4T6.5" "7.5 k4T7.5"; do
    set -- $spec
    bits=$1 name=$2
    [[ -s $OUT/index/$name.tidx ]] && continue
    echo "building $name"
    mkdir -p "$OUT/index"
    run filter build "$OUT/prep" 4 "$bits" "$OUT/index/$name.idx"
    run filter tile "$OUT/prep" "$OUT/index/$name.idx" "$OUT/index/$name.tidx.part"
    mv "$OUT/index/$name.tidx.part" "$OUT/index/$name.tidx"
    rm -f "$OUT/index/$name.idx"
done

echo "prep:      $OUT/prep"
echo "index A:   $OUT/index/k4T6.5.tidx"
echo "index B:   $OUT/index/k4T7.5.tidx"
echo "image:     $SIF"
