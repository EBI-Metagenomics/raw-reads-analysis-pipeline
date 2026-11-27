# Changelog

All notable changes to the Raw Reads Analysis Pipeline (RawR) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `pfam_orf_filter`: keep ORFs in FragGeneScanRs's frame, plus other ORFs with a hexamer coding score at or above `pfam_orf_filter_threshold`, before the Pfam search. Off by default.
- `pfam_search_engine = nail`: nail's seeding finds candidate (model, ORF) pairs, then HMMER checks only those pairs at gathering thresholds, through `pyhmmer`. Alignments are published as `*_pfam_alignments.parquet`. The default stays `hmmsearch`.
- `pfam_search_engine = kfpfam`: kfpfam's k-mer prefilter in place of nail's seeding, with HMMER checking its candidate pairs exactly as the nail path does. Its prep directory, tiled index and image are built once by `hpc/pfam-arms/kfpfam_setup.sh` and passed as `kfpfam_prep`, `kfpfam_index` and `kfpfam_sif`; `kfpfam_bar`, `kfpfam_gate`, `kfpfam_min_gapped` and `kfpfam_block` set the mode, defaulting to kfpfam's Mode A. The default engine stays `hmmsearch`.
- `pfam_seeded_chunk_reads` (4 M): splits the reads before translation for the `nail` and `kfpfam` engines, which search one whole HMM file and so have no database chunking to spread a sample over. Translation, seeding and HMMER are then one task per chunk rather than one per sample. Set it to 0 for one task per sample. The `hmmsearch` path is unchanged.
- `pfam_search_engine = nail` works with a chunked Pfam database as well as a single `Pfam-A.hmm`. nail searches one whole HMM file, so each chunk is seeded and searched separately, then the domtbls are concatenated and the alignments merged into one Parquet per sample. `-Z` stays the whole database's model count, so E-values match an unchunked search.

- `pfam_weak_evalue` (0, off): for the `nail` and `kfpfam` engines, also writes domains under Pfam's gathering thresholds to `*_pfam_alignments.parquet`, marked by the new `above_ga` column. The domtbl, the profile and its statistics don't change. On 12,000 windows of 50 residues and Pfam 38.2, E 0.0301 (P 1e-6) gives 3,556 domains more than the 4,115 above the thresholds, for 7 core-s of HMMER against 6.
- `pfam_lookup` (off): profiles each sample's Pfam alignments against an index of MGnify's proteins aligned to the same models, for the closest protein and its labels (`*_pfam_lookup.*`). `nail` and `kfpfam` only. The index (`pfam_lookup_bundle`) is built by kmer-functional-profiler's `workflows/lookup-index`, and the image (`pfam_lookup_sif`) by `hpc/pfam-arms/lookup_setup.sh`; alignments to other model versions than the index's are an error. `pfam_lookup_floor`, `pfam_lookup_k`, `pfam_lookup_second` and `pfam_lookup_margin` set the calls; `pfam_lookup_chunk_size` splits a sample's lookup into runs of families, merged by `PFAM_LOOKUP_MERGE`.
- `*_pfam_alignments.parquet` gains `pfam_version`, each model's version, so that a reader can refuse columns of another release.
- `pfam_alignments` (on): whether the `nail` and `kfpfam` engines write and publish each domain's HMM alignment. They go to zstd Parquet as they are found, about 55 bytes a domain, where they were held as Python objects (about 700 bytes) until the end of the task.

### Fixed

- `pfam_search_engine = kfpfam` with a prep directory built from another Pfam release is now an error. Seeds join to models by name, and names the HMM file lacks were skipped, so renamed families got no reads without a warning. The kfpfam task compares the number of models with `num_models`, `hmmer_pairs.py` refuses models that the HMM files lack, and `hpc/pfam-arms/kfpfam_setup.sh` rebuilds a prep made from a different file.
- `*_pfam_alignments.parquet` gains `read_from` and `read_to`, the aligned residues' nucleotides on the read, so a domain can be placed on its read without redoing the frame arithmetic. `hmmer_domtbl_parse_coverage.py`'s own mapping is not usable for this: it adds the frame offset a second time on top of the ORF's `begin` and places a reverse frame one nucleotide low, which cancels out when it compares hits on one read but is wrong as a read coordinate.

- The `hmmsearch_subsampling` draw is now reproducible. `SEQKIT_SHUFFLE_FASTA` picked its reads with unseeded `shuf`, so every run searched a different subsample and two runs of the same data could not be compared. It now uses seeded reservoir sampling (`ext.seed`, 23 by default), which also holds only the kept read names in memory rather than all of them.

### Changed

- `KFPFAM_FILTER` and the `HMMER_PAIRS` task after it are one task, `KFPFAM_HMMER`. The filter's candidate pairs go to HMMER through a pipe, a block of ORFs at a time, so no seeds are written (11 GB per 4 M reads at kfpfam's most sensitive setting, after as much again in a temporary file) and memory doesn't grow with the sample: `HMMER_PAIRS` held every candidate ORF, about 1 KB each, and peaked at 42 GB. Hits, scores and E-values are the same; the domain table's c-Evalue column, which nothing reads, depends on the batch. The kfpfam image must be rebuilt, as it now holds pyhmmer and pyarrow, and config selectors naming `KFPFAM_FILTER` must become `KFPFAM_HMMER`. `nail` still uses `HMMER_PAIRS`.
- `MERGEPARQUET` sorts a range of families at a time, from temporary files, rather than a whole sample: 1.9 GB at 36 M domains, where it took up to 27 GB.
- `subworkflows/local/nail_hmmer` is now `seeded_hmmer`, since nail and kfpfam share it: a seeding step, then HMMER on the candidate pairs. Config selectors naming `NAIL_HMMER` must become `SEEDED_HMMER`.
- Pfam moves from 38.0 (25,545 models) to 39.0 (53,275 models), so Pfam profiles change. Existing `pfam` database caches must be cleared or refreshed with `--force_download_dbs`: the cache is found by name, not by release, so a 38.0 or 38.2 one is otherwise searched with 39.0's `-Z`. A site config with its own copy or chunks of Pfam has to move with it.
- `bbmap/reformat` and `bbmap/reformat_standardise` replaced with the shared modules from [ebi-metagenomics/nf-modules](https://github.com/EBI-Metagenomics/nf-modules), with the channel and `modules.config` changes that required ([#26](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/26)).
- GitHub Actions now runs the nf-test suite on `master` under the Docker profile.
- End-to-end tests assert on `.fastp.json` content and line counts instead of md5 checksums, which were not stable across runs.

## [6.0.0] - 2026-04-13

First release of this pipeline, for MGnify v6. The MGnify v5 raw-reads analysis has been
rewritten in Nextflow DSL2: taxonomic profiling is a re-implementation of the v5 method,
functional profiling is new, and the surrounding configuration, testing and documentation
were built from scratch in this repository. Versioned to match the MGnify platform.

### Added

- Functional profiling against Pfam-A: `seqkit` translation of subsampled reads, `hmmsearch`, and coverage parsing of the HMMER domtbl output. Not part of MGnify v5.
- ENA-aligned samplesheet schema (`sample`, `fastq_1`, `fastq_2`, `single_end`, `instrument_platform`) with nf-schema parameter and samplesheet validation.
- Stage skip flags: `skip_qc`, `skip_decontam`, `skip_taxonomic` ([#23](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/23)), `skip_functional` and `skip_standardise`.
- Chunking and subsampling parameters for `cmsearch`, `hmmsearch` and decontamination, so long-read and large studies fit within job limits.
- Database handling: automatic download and caching, per-database local path overrides, `force_download_dbs`, and optional retrieval via the EBI `downloadfromfire` module.
- MultiQC reports per run and for the whole pipeline run.
- nf-test coverage for local modules, subworkflows and the whole pipeline, with miniature test databases and read sets, run in GitHub Actions ([#10](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/10)).
- `docker`, `singularity`, `conda`, `local` and `test` profiles, with resource limits, retry policies and resource escalation for the long-running search steps ([#2](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/2), [#3](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/3), [#4](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/4)).
- Documentation: README with parameter reference, pipeline schema diagram, `docs/output.md` and `CITATIONS.md`.

### Changed

- Whole pipeline re-implemented in Nextflow DSL2, replacing the MGnify v5 implementation, with modules and subworkflows taken from nf-core and [ebi-metagenomics/nf-modules](https://github.com/EBI-Metagenomics/nf-modules) wherever possible rather than written locally.
- Quality control now uses `fastp`, with read standardisation and paired-end repair by BBMap `reformat`/`repair`. Standardisation is on by default, as it is needed in practice with ENA data.
- Decontamination (phiX and host hg38) re-implemented following ebi-metagenomics/nf-modules: `bwa-mem2` for short reads, `minimap2` for long reads, chunked, with alignment and filtering merged so BAM files are never stored.
- Taxonomic profiling re-implemented: mOTUs, and SILVA SSU/LSU via Rfam `cmsearch` followed by MAPseq, `mapseq2biom` and Krona plots.
- Long reads (`OXFORD_NANOPORE`, `PACBIO_SMRT`) and single-, paired-end and interleaved short reads are each routed through QC, decontamination and standardisation differently, instead of the single short-read path of v5.
- Outputs renamed, split into `qc_passed` and `qc_failed` summaries, gzipped and given headers, for orchestration by Prefect within the MGnify system.
- Read staging handled by Nextflow rather than custom fetch modules.
- Per-run MultiQC moved to a local module that takes `meta` in its input, and the MultiQC module updated from 1.25 to 1.33 ([#19](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/19), [#17](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/17)).
- Codon-specific profile removed, its settings becoming pipeline defaults.

### Fixed

Issues found and fixed during development of this release, before it was published.

- Long-read instrument platform names now match ENA, and more short-read platforms are supported ([#5](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/5)).
- Phred-33 quality encoding is now assumed explicitly; BBMap auto-detection misread Oxford Nanopore data as Phred-64 and failed ([#6](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/6)).
- Invalid paired-end FASTQ files handled by adding `bbmap_repair` to drop unpaired reads ([#18](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/18)).
- `hmmsearch` and MultiQC failures under Docker, including on macOS ([#20](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/20)).
- HMMER output parsing: zero-length matches are no longer produced, and query coordinates are corrected for the off-by-one that made length-1 matches appear empty ([#24](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/24)).
- Repeated headers in the HMM output ([#21](https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/pull/21)).
- Sequence coverage calculation for Pfam mapping, subsampling that took the head of the file rather than shuffling, and empty SILVA SSU/LSU results stopping the run.
- Output filenames: trailing dots, missing `.gz` suffixes, and inputs being overwritten by outputs during decontamination.

[Unreleased]: https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/compare/v6.0.0...master
[6.0.0]: https://github.com/EBI-Metagenomics/raw-reads-analysis-pipeline/releases/tag/v6.0.0
