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
[[ -f run.sh && -f site.config ]] || { echo "submit from hpc/pfam-arms" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g' NXF_VER=25.10.4
PIPELINE=../..

arm() {
    local name=$1
    shift
    nextflow run "$PIPELINE" -profile singularity -c site.config -resume \
        --samplesheet samplesheet.csv --outdir "results/$name" \
        --skip_functional false --skip_taxonomic true \
        --trace_report_suffix arm --singularity_cachedir singularity_cache "$@"
}

arm current
arm filter --pfam_orf_filter true
arm nail --pfam_search_engine nail
arm both --pfam_search_engine nail --pfam_orf_filter true

: > results/compare.tsv
: > results/families.tsv
for sample in $(tail -n +2 samplesheet.csv | cut -d, -f1); do
    args=()
    for a in current filter nail both; do
        args+=("$a=results/$a/$sample/function-summary/pfam/${sample}_pfam.txt.gz")
    done
    python3 "$PIPELINE/bin/compare_profiles.py" "$sample" "${args[@]}" --families "results/$sample.families.tsv" >> results/compare.tsv
    cat "results/$sample.families.tsv" >> results/families.tsv
    rm "results/$sample.families.tsv"
done
python3 cost.py current=results/current/pipeline_info/execution_trace_arm.txt filter=results/filter/pipeline_info/execution_trace_arm.txt \
    nail=results/nail/pipeline_info/execution_trace_arm.txt both=results/both/pipeline_info/execution_trace_arm.txt > results/cost.tsv
echo finished
