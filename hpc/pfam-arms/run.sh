#!/usr/bin/env bash
#SBATCH --job-name=pfam-arms
#SBATCH --output=%x.%j.log
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# V5 of docs/plans/pfam-orf-filter-nail.md: four arms of Pfam profiling on samplesheet.csv.
# Submit from this directory: cd hpc/pfam-arms && sbatch run.sh
# The standard run, kfpfam top with the MGnify lookup (lookup_setup.sh, then kfpfam_setup.sh):
#   RUN_ARMS=kfpfam_top sbatch run.sh
# CONFIGS=<pipelines-nextflow-configs checkout>: Codon's config for this pipeline too (its local host, phiX,
# Rfam and SILVA databases, so nothing is downloaded), under site.config; the work directory stays work/ here.
# The reads come from ENA's FIRE (--use_fire_download, as B11): samplesheet.csv's paths have no scheme, so
# without it they are read as local files. It needs the Nextflow secrets FIRE_ACCESS_KEY and FIRE_SECRET_KEY
# (`nextflow secrets list`); FIRE_DOWNLOAD=false for a samplesheet of local or ftp:// paths.
# Re-submitting resumes. Each arm reuses the others' QC, decontamination and translation.
# Read:
#   results/compare.tsv   per sample and arm, against `current`: read and family differences, worst families
#   results/compare_weefams.tsv  the same for kfpfam_top_weefams against current_weefams (weefams_setup.sh)
#   results/families.tsv  reads per family in each arm
#   results/cost.tsv      core-hours and peak memory per Pfam step
#   results/<arm>/<sample>/function-summary/pfam/  each arm's published files, with filter stats
set -euo pipefail
[[ -f run.sh && -f site.config && -f unchunked_pfam.config ]] || { echo "submit from hpc/pfam-arms" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g' NXF_VER=25.10.4
PIPELINE=${PIPELINE:-../..}  # a pipeline checkout; this one by default

# The nail arms search one whole Pfam-A.hmm, where the hmmsearch arms use the site's chunks.
# Set PFAM_HMM to a single Pfam-A.hmm if you have one; otherwise set PFAM_CHUNKS to the chunk
# files the site config serves (a glob is fine) and this builds it once, beside this script.
# An HMM file is its models one after another, so concatenating the chunks gives the database
# the hmmsearch arms search, and -Z below is counted from the result rather than assumed.
# Only the nail arms need it: RUN_ARMS without them (kfpfam_top alone) skips this.
nail_arms=1
if [[ -n ${RUN_ARMS:-} && " $RUN_ARMS " != *" nail "* && " $RUN_ARMS " != *" nail_loose "* && " $RUN_ARMS " != *" both "* ]]; then
    nail_arms=0
fi
PFAM_HMM=${PFAM_HMM:-$PWD/pfam_unchunked/Pfam-A.hmm}
if (( nail_arms )) && [[ ! -s $PFAM_HMM ]]; then
    [[ -n ${PFAM_CHUNKS:-} ]] || {
        echo "set PFAM_HMM to one Pfam-A.hmm, or PFAM_CHUNKS to the chunk HMMs to join" >&2
        exit 1
    }
    mkdir -p "$(dirname "$PFAM_HMM")"
    cat $PFAM_CHUNKS > "$PFAM_HMM.part"
    mv "$PFAM_HMM.part" "$PFAM_HMM"
fi
if (( nail_arms )); then
    export PFAM_HMM_DIR=$(dirname "$PFAM_HMM")
    export PFAM_HMM_NAME=$(basename "$PFAM_HMM")
    export PFAM_NUM_MODELS=$(grep -c '^NAME ' "$PFAM_HMM")
    echo "nail arms search $PFAM_HMM, $PFAM_NUM_MODELS models"
fi

codon=()
if [[ -n ${CONFIGS:-} ]]; then
    codon=(-c "$CONFIGS/pipelines/raw-reads-analysis-pipeline/codon.config")
    [[ -s ${codon[1]} ]] || { echo "no ${codon[1]}: CONFIGS is not a pipelines-nextflow-configs checkout" >&2; exit 1; }
fi

# RUN_ARMS: the arms to run, by name ('kfpfam_top' for the standard run alone); all of them if unset
arm() {
    local name=$1
    shift
    [[ -z ${RUN_ARMS:-} || " $RUN_ARMS " == *" $name "* ]] || return 0
    nextflow run "$PIPELINE" -profile singularity ${codon[@]+"${codon[@]}"} -c site.config -resume -work-dir work \
        --samplesheet samplesheet.csv --outdir "results/$name" \
        --skip_functional false --skip_taxonomic true \
        --hmmsearch_subsampling=-1 \
        --trace_report_suffix arm --singularity_cachedir singularity_cache \
        --use_fire_download "${FIRE_DOWNLOAD:-true}" "$@"
}

# Seeded arms first: they are the cheap ones, and whichever runs first pays for the shared QC,
# decontamination and merge, so the hmmsearch arms start from a warm cache rather than the other
# way round. With subsampling off, decide from the seeded arms' real cost whether the hmmsearch
# arms are worth running at full depth at all. kfpfam goes first of all: its prefilter is the
# cheapest, so a failure there surfaces before anything expensive is queued.

# kfpfam's top mode, the standard for reads, with the MGnify lookup on its alignments: on the Pfam its
# index was built from (kfpfam/Pfam-A.hmm, kfpfam_setup.sh), which is the lookup bundle's own, Pfam 38.1.
# The pipeline's kfpfam defaults are top's (bar 15, gate 14.5, HMMER at 21.25), and KFPFAM_HMMER takes the
# gapped bars from the prep. One whole Pfam, like the nail arms: the index covers the whole database. Read
# chunks of 10 M, not the pipeline's 4 M: a 4 M-read task took about 15 minutes, and memory stays flat.
if [[ -s kfpfam/index/top.tidx ]]; then
    kf="--kfpfam_sif $PWD/kfpfam/kfpfam-filter.sif --kfpfam_prep $PWD/kfpfam/prep --kfpfam_index $PWD/kfpfam/index/top.tidx"
    # The bundle must be of the Pfam searched here, or PFAM_LOOKUP stops. Everything before the search is cached.
    if [[ -s lookup/bundle/manifest.json && -s lookup/kfp-lookup.sif ]]; then
        kf+=" --pfam_lookup true --pfam_lookup_bundle $(realpath lookup/bundle) --pfam_lookup_sif $PWD/lookup/kfp-lookup.sif"
    else
        echo "no lookup/bundle and lookup/kfp-lookup.sif: the kfpfam arm runs without the MGnify lookup" >&2
    fi
    PFAM_HMM_DIR=$PWD/kfpfam PFAM_HMM_NAME=Pfam-A.hmm PFAM_NUM_MODELS=$(grep -c '^NAME ' kfpfam/Pfam-A.hmm) \
        arm kfpfam_top -c unchunked_pfam.config $kf --pfam_search_engine kfpfam --pfam_seeded_chunk_reads 10000000
else
    echo "no kfpfam/index/top.tidx: run kfpfam_setup.sh to add the kfpfam arm" >&2
fi

# Weefams (weefams_setup.sh): kfpfam's top mode on reads, and hmmsearch on every ORF as its truth. hmmsearch costs
# what the library's columns cost against Pfam's (39.0 holds about 10 M), so look at the weefams' count first.
if [[ -s weefams/index/top.tidx ]]; then
    export WEEFAMS_DIR=$PWD/weefams WEEFAMS_MODELS=$(grep -c '^NAME ' weefams/weefams.hmm)
    echo "weefam arms search $WEEFAMS_MODELS weefams"
    arm kfpfam_top_weefams -c weefams.config --kfpfam_sif "$PWD/kfpfam/kfpfam-filter.sif" --pfam_search_engine kfpfam
    arm current_weefams -c weefams.config
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

[[ -s results/current/pipeline_info/execution_trace_arm.txt ]] || { echo "no current arm: nothing to compare against"; echo finished; exit 0; }
# `current` first, so compare_profiles.py measures every other arm against it
ARMS="current filter nail both nail_loose"
if [[ -s kfpfam/index/top.tidx ]]; then
    ARMS="$ARMS kfpfam_top"
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
# the weefam arms against each other: weefams are not Pfam's families
WEE=
[[ -s weefams/index/top.tidx ]] && WEE="current_weefams kfpfam_top_weefams"
if [[ -n $WEE ]]; then
    : > results/compare_weefams.tsv
    for sample in $(tail -n +2 samplesheet.csv | cut -d, -f1); do
        args=()
        for a in $WEE; do args+=("$a=results/$a/$sample/function-summary/pfam/${sample}_pfam.txt.gz"); done
        python3 "$PIPELINE/bin/compare_profiles.py" "$sample" "${args[@]}" >> results/compare_weefams.tsv
    done
fi
cost=()
for a in $ARMS $WEE; do cost+=("$a=results/$a/pipeline_info/execution_trace_arm.txt"); done
python3 cost.py "${cost[@]}" > results/cost.tsv
echo finished
