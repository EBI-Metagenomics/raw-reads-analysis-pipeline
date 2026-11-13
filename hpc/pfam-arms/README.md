# Pfam arms: V5 of the ORF filter and nail plan

Runs the pipeline four times on `samplesheet.csv`: `current` (`hmmsearch` on all ORFs), `filter`, `nail` and `both`. Parameters other than these two options stay at the pipeline's defaults, including the cap of 2 M reads for Pfam, except that taxonomic profiling is skipped.

```bash
cd hpc/pfam-arms && PFAM_CHUNKS='/path/to/pfam/chunks/*/Pfam-A.hmm' sbatch run.sh
```

The nail arms search one whole `Pfam-A.hmm`, not the site's chunks, so `run.sh` needs to know where
to get one. Set `PFAM_HMM` if a single file already exists, or `PFAM_CHUNKS` to the site's chunk
files and it joins them once into `pfam_unchunked/Pfam-A.hmm` beside the script; an HMM file is just
its models one after another, so that is the same database the hmmsearch arms search. `-Z` is
counted from the result rather than assumed, so the arms stay comparable. Neither variable is needed
on a re-submission once the file exists.

The first run downloads the host and phiX databases, and the FASTQ files from ENA. On one local test, `hmmsearch` used about 60 core-ms per read, so `current` could be about 32 core-hours per sample at the 2 M cap. Local timings have run well below the HPC's, so allow several times that.

## Chunking and resources

`site.config` sizes the search steps from B11's traces, where a task searched one sample's whole ORF
set against one Pfam chunk:

| step | tasks per sample | median | worst | cores used |
|---|---|---|---|---|
| `HMMER_HMMSEARCH` | 20 | 65 min | 144 min | 2.25 of 4 |
| `NAIL_SEED` | 20 | 1.5 min | 3.9 min | 3.7 of 16 |
| `HMMER_PAIRS` | 20 | 36 s | 1.3 min | 2.1 of 8 |

So the two engines wanted opposite changes. `hmmsearch` tasks were long and used barely two cores,
and the reads were never split, since the default `hmmsearch_chunksize` of 1 GB exceeds a sample's
ORFs: `params.hmmsearch_chunksize = 256.MB` makes four read chunks of about a quarter the length
each, and `cpus = 2` stops the other two sitting idle. nail's tasks were slivers, each paying to
rebuild its MMseqs2 query database over the whole ORF set, so the nail arms now search an unchunked
Pfam and run one `NAIL_SEED` and one `HMMER_PAIRS` per sample.

Memory for the unchunked `NAIL_SEED` is the one figure B11 cannot give: 16.6 GB was its peak against
a twentieth of Pfam, and the database side is what grows, so `site.config` asks for 48 GB and
escalates on retry. Check `results/cost.tsv` afterwards and bring it down to what the run used. There
is no `disk` request: nail's `--tmp-dir` is the task's work directory on `/hps/nobackup`, so asking
Slurm for node-local space only stops the job being scheduled.

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
