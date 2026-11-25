#!/usr/bin/env bash
#SBATCH --job-name=lookup-setup
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=lookup-setup.%j.log
#
# One-off setup for the kfpfam arms' MGnify lookup: the image, and a link to the bundle. Run it on the
# HPC from this directory, before run.sh. Everything lands in lookup/ beside this script.
#
#   BUNDLE=<lookup-index results>/bundle ./lookup_setup.sh
#
#   BUNDLE  the bundle of kmer-functional-profiler's workflows/lookup-index, for the Pfam the arms search
set -euo pipefail
: "${BUNDLE:?the bundle directory of lookup-index}"
[[ -s $BUNDLE/manifest.json ]] || { echo "$BUNDLE/manifest.json is not there" >&2; exit 1; }
OUT=$PWD/lookup
SIF=$OUT/kfp-lookup.sif
mkdir -p "$OUT"

# Published nowhere: built where Docker is, then copied here as a tar, as kfpfam's image is. Rebuilt when
# the tar is newer, so a re-exported image doesn't leave a stale .sif behind.
if [[ ! -s $SIF || ( -s $OUT/kfp-lookup.tar && $OUT/kfp-lookup.tar -nt $SIF ) ]]; then
    if [[ ! -s $OUT/kfp-lookup.tar ]]; then
        cat >&2 <<MSG
no image. Build it where Docker is, from a kmer-functional-profiler checkout, and copy the tar to $OUT:
  docker build --platform linux/amd64 -t kfp-lookup workflows/pfam-anchor/lookup
  docker save kfp-lookup:latest -o kfp-lookup.tar
Then re-run this script.
MSG
        exit 1
    fi
    rm -f "$SIF"
    singularity build "$SIF" "docker-archive://$OUT/kfp-lookup.tar"
fi

ln -sfn "$(realpath "$BUNDLE")" "$OUT/bundle"
singularity exec --bind "$(realpath "$BUNDLE")" "$SIF" python3 -c "
import json, sys
m = json.load(open(sys.argv[1] + '/manifest.json'))
print('bundle: format', m['format'], 'Pfam', m['pfam'], 'MGnify', m['mgnify'], m['rows'], 'rows')
" "$(realpath "$BUNDLE")"
echo "image:  $SIF"
