# Changelog

All notable changes to the Raw Reads Analysis Pipeline (RawR) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

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
