# Pfam arms: V5 of the ORF filter and nail plan

Runs the pipeline on `samplesheet.csv` once per arm: `current` (`hmmsearch` on all ORFs), `filter`, `nail` and `both`. Parameters other than these options stay at the pipeline's defaults, except that taxonomic profiling is skipped and Pfam subsampling is off (`--hmmsearch_subsampling -1`), so every merged read is searched.

```bash
cd hpc/pfam-arms && PFAM_CHUNKS='/path/to/pfam/chunks/*/Pfam-A.hmm' sbatch run.sh
```

The nail arms search one whole `Pfam-A.hmm`, not the site's chunks, so `run.sh` needs to know where
to get one. Set `PFAM_HMM` if a single file already exists, or `PFAM_CHUNKS` to the site's chunk
files and it joins them once into `pfam_unchunked/Pfam-A.hmm` beside the script; an HMM file is just
its models one after another, so that is the same database the hmmsearch arms search. `-Z` is
counted from the result rather than assumed, so the arms stay comparable. Neither variable is needed
on a re-submission once the file exists.

The first run downloads the host and phiX databases (unless `CONFIGS` is set), and the FASTQ files from ENA's FIRE (`--use_fire_download`, which needs the Nextflow secrets `FIRE_ACCESS_KEY` and `FIRE_SECRET_KEY`; the samplesheet's paths have no scheme, so they are not read without it). B11 measured, at the 2 M cap over eleven samples, 550 core-hours for `current`, 358 for `filter`, 27 for `nail` and 17 for `both`. Subsampling is off now and the eleven libraries hold about twenty times those reads, so scale accordingly: `current` alone is of the order of 10,000 core-hours. The arms run seeded-first for that reason, kfpfam ahead of nail — read the cheap arms' real numbers before committing to the hmmsearch pair at full depth.

## Chunking and resources

`site.config` sizes the search steps from B11's traces, where a task searched a 2 M-read sample's
whole ORF set against one Pfam chunk. Those numbers are a floor now: subsampling is off, so a
library holds about twenty times the reads.

The seeded engines search one whole HMM file, so they have no database chunking to spread a sample
over. `pfam_seeded_chunk_reads` (4 M) splits the reads instead, before translation, so translation,
seeding and HMMER are each one task per 4 M reads rather than one per sample: about ten tasks per
library here. Set it to 0 to go back to one task per sample. The hmmsearch arms are untouched —
they keep chunking ORFs by size, since they already spread over the database's chunks.

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

## The standard run: kfpfam top and the MGnify lookup

kfpfam's top mode is the standard for reads, and its alignments go straight to the MGnify lookup
(`--pfam_lookup`). It searches the Pfam the lookup bundle was built from, Pfam 38.1, so the reads'
alignments and the bundle agree; the other arms keep the site's Pfam. Once per bundle, from this directory:

```bash
mkdir -p lookup kfpfam
cp <kfp runs>/kfp-ablations/images/kfp-lookup.tar lookup/
BUNDLE=<kfp runs>/kfp-ablations/lookup-index/results/bundle sbatch lookup_setup.sh
cp <kfp runs>/kfp-ablations/images/kfpfam-filter.tar kfpfam/
KFPFAM=<kfpfam, prep-any-hmm> KFP=<kmer-functional-profiler, pfam-anchor> sbatch kfpfam_setup.sh
CONFIGS=<pipelines-nextflow-configs checkout> RUN_ARMS=kfpfam_top sbatch run.sh
```

`kfpfam_setup.sh` builds top from `lookup/bundle/hmm/Pfam-A.hmm` (or `PFAM_HMM`) into `kfpfam/`, with a copy
of that file, which the arm searches. A prep of another file (39.0's, for the old A and B arms) is moved
to `kfpfam/stale.<date>/`. Top's settings are the pipeline's kfpfam defaults: seeds of shape 10110001 at
theta 9, bar 15, gate 14.5, HMMER at 21.25 or the model's gapped bar (`prep/gapped.txt`, which
`KFPFAM_HMMER` uses when it is there), band 15 at any diagonal. The arm publishes
`results/kfpfam_top/<sample>/function-summary/pfam/`: `*_pfam_alignments.parquet` and
`*_pfam_lookup.{reads.parquet,ko.tsv,summary.json}`. `CONFIGS` adds Codon's config for this pipeline under `site.config`, so the host, phiX, Rfam and SILVA
databases are its local copies rather than downloads; the work directory stays `work/` here. Only the nail arms
need `pfam_unchunked/Pfam-A.hmm` (or `PFAM_CHUNKS`). `RUN_ARMS` names the arms to run (all if unset); the
comparison at the end needs `current`, and is skipped without it. Against `current`, kfpfam_top is a
different Pfam release (38.1 against the site's), so read its `compare.tsv` rows with that in mind.

`PFAM_LOOKUP` stops at the first family of another version than the bundle's. To profile alignments that
exist already without searching again, use kmer-functional-profiler's `hpc/kfp-ablations/lookup-profile`.

## Samples

`samplesheet.csv` and `runs.tsv` are the response's own reads, copied from
`nar-review-response/hpc/B11_pfam_arms`, so the Zymo runs have ground truth and the biome runs are
comparable with B2. Six Zymo standards cover the mock communities, five natural metagenomes the
biomes, and two further runs sit at 250 nt for the read-length question; the `long` pair is run
separately, through `samplesheet_long.csv` there.

| Run | Set | Label | Reads | Read length |
|---|---|---|---|---|
| ERR2984773 | zymo | D6300 even (cells) | 8.8 M | 151 |
| SRR12324253 | zymo | D6300 even (cells) | 51.5 M | 150 |
| SRR12829160 | zymo | D6305/D6306 even (DNA) | 40.7 M | 151 |
| ERR2935805 | zymo | D6310 log (cells) | 47.8 M | 101 |
| SRR23926920 | zymo | D6331 gut | 37.7 M | 150 |
| SRR32629731 | zymo | D6331 gut | 10.7 M | 99 |
| ERR6785304 | biome | bird gut | 73.5 M | 72 |
| ERR4334731 | biome | freshwater | 62.9 M | 145 |
| ERR5006025 | biome | human gut | 42.6 M | 112 |
| SRR12479991 | biome | marine | 23.5 M | 151 |
| ERR5866877 | biome | soil | 40.2 M | 75 |
| DRR160704 | long | soil | 9.9 M | 250 |
| DRR092732 | long | deep sea | 6.6 M | 251 |

Read length matters here because the ORF filter's coding score sums over an ORF, so a threshold set
on 150 nt reads need not hold at 72 or 250. `compare.tsv` is where to look for that.

## What to look at

- `results/compare.tsv`: whether `of_those_below_0.8` and `worst` hold up across samples, and whether they get worse with read length, because coding scores add up over an ORF.
- `results/cost.tsv`: core-hours and peak memory per step, to set `conf/modules.config`'s resources for `NAIL_SEED` and `HMMER_PAIRS`.
- `results/<arm>/<sample>/function-summary/pfam/*_pfam_orf_filter.stats.json`: share of residues kept, and whether a sample fell back to all ORFs.

## Weefams

Two more arms profile the reads with mgnifams-explore's weefams instead of Pfam, once `weefams/` exists:
`kfpfam_top_weefams` (kfpfam's top mode, as on reads) and `current_weefams` (`hmmsearch` on every ORF, its truth).

```bash
WEEFAMS=<runs>/weefams/full/results/wcollect/weefams.hmm KFP=<kmer-functional-profiler checkout> sbatch weefams_setup.sh
```

It needs `kfpfam/kfpfam-filter.sif` from `kfpfam_setup.sh`, of kfpfam 75e543d or later. It writes `weefams/weefams.hmm` with
an `ACC` line for each weefam (its name: the pipeline counts reads by accession, and weefams have none), the prep with
top's gapped bars, and `weefams/index/top.tidx`, built as hm-annotate's PREP builds top. `weefams.config` points the
Pfam step at them; the MGnify lookup is off, as its bundle is of Pfam. `results/compare_weefams.tsv` compares the two arms.
`current_weefams` costs what the library's columns cost against Pfam's, so check the number of weefams before it runs.
