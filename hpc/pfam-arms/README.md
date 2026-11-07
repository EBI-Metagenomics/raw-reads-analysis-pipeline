# Pfam arms: V5 of the ORF filter and nail plan

Runs the pipeline four times on `samplesheet.csv`: `current` (`hmmsearch` on all ORFs), `filter`, `nail` and `both`. Parameters other than these two options stay at the pipeline's defaults, including the cap of 2 M reads for Pfam, except that taxonomic profiling is skipped.

```bash
cd hpc/pfam-arms && sbatch run.sh
```

The first run downloads Pfam 38.2, the host and phiX databases, and the FASTQ files from ENA. On one local test, `hmmsearch` used about 60 core-ms per read, so `current` could be about 32 core-hours per sample at the 2 M cap. Local timings have run well below the HPC's, so allow several times that.

## Samples

Chosen from ENA's metadata only: paired Illumina, labelled metagenomic WGS, 5–25 M read pairs. Contents haven't been checked.

| Sample | Environment | Read length | Note |
|---|---|---|---|
| DRR012573 | Soil | 100 |  |
| DRR1002327 | Soil | 150 |  |
| DRR160704 | Soil | 250 |  |
| DRR083064 | Marine plankton | 100 | 20–100 µm filter fraction, so likely eukaryote-rich |
| SRR23559199 | Marine coastal | 150 |  |
| DRR092732 | Deep sea | 250 | 4,500 m |
| SRR30690279 | Freshwater | 150 |  |
| SRR28214405 | Viral metagenome | 150 | ENA's library source says GENOMIC |
| DRR1003640 | Human gut | 150 | Host reads removed by the pipeline |

The set doesn't deliberately cover high and low GC, because ENA has no GC field. `compare.tsv` shows where each sample lands; check the extremes among these before relying on the threshold.

## What to look at

- `results/compare.tsv`: whether `of_those_below_0.8` and `worst` hold up across samples, and whether they get worse with read length, because coding scores add up over an ORF.
- `results/cost.tsv`: core-hours and peak memory per step, to set `conf/modules.config`'s resources for `NAIL_SEED` and `HMMER_PAIRS`.
- `results/<arm>/<sample>/function-summary/pfam/*_pfam_orf_filter.stats.json`: share of residues kept, and whether a sample fell back to all ORFs.
