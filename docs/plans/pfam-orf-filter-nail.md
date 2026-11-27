# Plan: ORF filter and nail + HMMER for Pfam profiling

Branch `pfam-orf-filter-nail`, from `origin/master`. Both features are off by default. Pfam moves to 38.2 for every run, so default output changes only through that commit.

## Goal

1. **Pfam 38.2 everywhere**, the release behind MGnify's protein annotations (`mgy_proteins_pfam`).
2. **ORF filter:** keep ORFs in FragGeneScanRs's (FGS's) frame, plus other ORFs with hexamer coding score ≥ −5.
3. **Search engine `nail`:** nail finds candidate (model, ORF) pairs. `pyhmmer` then checks only those pairs at `--cut_ga` thresholds, and its alignments are stored in Parquet.
4. **Compare nail's and HMMER's alignments during development.** If they are the same, drop the HMMER step and use nail alone.

```
pfam_reads ─┬─ SEQKIT_TRANSLATE ───┐
            └─ FRAGGENESCANRS ─────┴─ SELECTFRAMES ─┐      (filter on; otherwise ORFs pass straight through)
                                                    │
      engine = hmmsearch:  splitFasta ─ HMMER_HMMSEARCH ─ COMBINEHMMSEARCHTBL ─┐
      engine = nail:       NAIL_SEED ─ HMMER_PAIRS ────────────────────────────┤
                                     └─ alignments.parquet                     │
                                                       PARSEHMMSEARCHCOVERAGE ─┴─ profile
```

`PARSEHMMSEARCHCOVERAGE` does not change. It reads a domtbl and needs seqkit's `_frame=_begin=_end=` names, so every new step keeps those names.

## Evidence so far

All on one sample, `ERR10889147`: 50,000 unmerged reads (309,149 ORFs), on a laptop at up to 6 threads, with Pfam 38.2 (30,134 models). Recall is against `hmmsearch` at Z = 20 M, E ≤ 1000, binned by best-domain score relative to GA.

**Candidates, with `pyhmmer` reporting down to E ≤ 1000:**

| Candidates | Pairs per read | ≥ GA | GA−5 to GA | −10 to −5 | −20 to −10 | < −20 | nail core-s | pyhmmer core-s |
|---|---|---|---|---|---|---|---|---|
| nail full run, `-S 1e-2`, reported hits | 10.7 | .989 | .908 | .832 | .760 | .681 | 701 | 15 |
| nail `--only-seed -S 1e-2`, MMseqs2 15 | 64.5 | .989 | .908 | .832 | .760 | .681 | not timed alone | 26 |
| nail `--only-seed -S 1e-2`, MMseqs2 18 | 65.8 | **.990** | **.909** | **.836** | **.768** | **.685** | 400 | 27 |
| nail `--only-seed -S 1e-1`, MMseqs2 15 | 334 | .990 | .909 | .833 | .762 | .686 | 548 | 65 |

Nail's times include its fixed cost of about 2 minutes wall. For scale, nail at default settings reached .987 at GA and `hmmsearch` costs 40–70 core-ms per read. Cost at scale is not measured yet (V4).

**`pyhmmer` at gathering thresholds** on the `-S 1e-2` seeds, `Z` = 30,134:

- It reported 8,447 domains on 8,366 (ORF, model) pairs in 11 s wall.
- Every pair is in the `hmmsearch` table, and all have the same best-domain score within 0.1 bits.
- Recall in the ≥ GA bin is .988 against .989 without cutoffs, and 37 pairs pass from below GA1. Both effects are expected from HMMER's separate sequence and domain thresholds, and V1 checks this.
- Parquet takes about 57 bytes per domain. At 0.17 domains per read that is about 19 MB per 2 M reads.

**nail's alignments against `hmmsearch`'s best domain** (preview from existing tables, nail full run at `-S 1e-2`). Coordinates are compared; alignment columns are not yet:

| Bin | Pairs | Same coordinates | Within 3 positions | nail − HMMER score, median (5–95%) |
|---|---|---|---|---|
| ≥ GA | 8,342 | 98.6% | 100% | 0.0 (−0.1 to 0.5) |
| GA−5 to GA | 3,070 | 94.6% | 99.8% | 0.1 (−0.1 to 0.7) |
| −10 to −5 | 7,304 | 89.7% | 99.8% | 0.2 (−0.2 to 0.8) |
| −20 to −10 | 21,176 | 86.9% | 99.6% | 0.3 (−0.3 to 1.0) |

Calling GA on best-domain score, the two agree on 8,328 pairs, with 78 passing only on nail's score and 103 only on HMMER's. So they are close but not the same, about 2% of GA calls differ. Nail alone is also not cheaper at this setting: 701 core-s against about 430 for seeds plus `pyhmmer`, though those two runs used different MMseqs2 versions.

What this settles:

- **Checking only the candidate pairs works in `pyhmmer`,** including at gathering thresholds.
- **Use seeds as candidates.** Recall matches nail's reported hits, at lower cost.
- **`-S 1e-2` is the permissive default.** `-S 1e-1` gives five times the pairs and no clear gain. The remaining 1% at GA comes from seeding itself.
- **MMseqs2 version:** the nail biocontainer has both `nail` 0.7.1 and MMseqs2 18.8cc5c, whereas nail is tested against 15-6f452. Version 18 did at least as well here, tested with the same MMseqs2 version in `kfpfam-tools`.

Not settled:

- **Exact match to `hmmsearch --cut_ga`** with the pipeline's flags (V1).
- **Alignment columns** between nail and HMMER (V2).
- **The filter and nail together:** their losses may add up (V3).
- **Filter inputs:** the filter was tested on unmerged reads with the `kfpfam` ORFs, not merged fragments with seqkit ORFs (V3).

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `pfam_orf_filter` | `false` | Run FGS and `SELECTFRAMES` |
| `pfam_orf_filter_threshold` | `-5` | Coding score for ORFs outside FGS's frame |
| `pfam_search_engine` | `hmmsearch` | `hmmsearch` or `nail` |
| `pfam_nail_seed_pvalue` | `1e-2` | nail `-S` |

The hexamer training size is fixed in the script, not a parameter, because −5 is only calibrated for that size.

## Steps

### 1. Pfam 38.2

Its own commit, so it can be reviewed or reverted alone. In `nextflow.config`, under `databases.pfam`:

```
remote_path    = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz"
remote_checksum = "7ab3c4e215d0daaea3004e37c4e24f8a"   // from the release's md5_checksums; replaces remote_checksum_path
files { hmm = "Pfam-A.hmm" }
variables { num_models = 30134 }
```

- `FETCHUNZIP` unpacks this to `pfam/Pfam-A.hmm`. An existing `pfam` cache directory has only the 38.0 file, so runs must use `--force_download_dbs` or clear it.
- Update the table in `README.md` and add a `CHANGELOG.md` entry stating that default profiles change.
- `conf/test.config` keeps `pfam_mini`, so its snapshots are unaffected. It sets `remote_checksum = false`, because it would otherwise inherit 38.2's checksum.

### 2. `HMMER_PAIRS`

`modules/local/hmmer_pairs/` with `bin/hmmer_pairs.py`, using `pyhmmer` 0.12.3 and `pyarrow` in a Seqera Community container, as `minimap2_align_samtools_bam2fq` does.

- **Inputs:** ORF FASTA, nail seeds, Pfam HMM, `num_models`.
- **Search:** read the HMM file once. For each model with candidates, build a sequence block from that model's targets and run `Pipeline(alphabet, Z=num_models, bit_cutoffs="gathering").search_hmm(...)` with default filters, on threads from `task.cpus`.
- **Output selection:** write only hits and domains marked `reported`, as `hmmsearch` does in its domtbl.
- **Memory:** hold candidates as one integer array per model and load only ORFs that have candidates. A Python set per model, as in the test, won't fit 2 M reads at about 130 M pairs. One task per sample. `# ponytail:` split by ORF chunk if V4 shows memory is a problem.
- **Model names:** seeds name the model by NAME even with `--use-accession`, so match on `hmm.name`.
- **Output:** `*.domtbl.gz` in `hmmsearch` column order and `*_alignments.parquet`.
- **Empty input:** no candidates gives an empty domtbl and a Parquet file with only the schema.

Parquet, zstd level 9, sorted by `pfam` then `read`, one row per GA-passing domain:

| Columns | Type | Notes |
|---|---|---|
| `read`, `frame`, `orf_begin`, `orf_end` | string, int8, int32, int32 | Parsed from the ORF name |
| `pfam`, `domain` | dictionary string, int16 | Accession without version |
| `seq_score`, `dom_score` | float32 | Bits |
| `seq_evalue`, `i_evalue` | float64 | At `Z = num_models` |
| `hmm_from`, `hmm_to`, `ali_from`, `ali_to`, `env_from`, `env_to` | int32 | |
| `cigar` | string | M, I, D against the HMM |
| `residues` | string | Aligned ORF residues, ungapped |
| `pp` | string | HMMER's posterior characters per aligned column |

Model consensus is left out: it can be rebuilt from `pfam`, `hmm_from` and `cigar` given Pfam 38.2. A check confirms the round trip from `cigar` and `residues` to HMMER's aligned sequence. Alignments below GA are not stored, because `--cut_ga` drops them.

### 3. `NAIL_SEED`

`modules/local/nail_seed/`, stock `nail:0.7.1--hab7d0fd_0`:

```
nail search --only-seed -S ${params.pfam_nail_seed_pvalue} -t ${task.cpus} --tmp-dir tmp --seeds-out ${prefix}.seeds.tsv ${hmm} ${orfs}
```

It runs once per sample, because each run has a fixed cost of about 2 minutes. Delete `tmp` in the script. Expect about 5 GB of uncompressed seeds per 2 M reads, written to the work directory only, unless the filter reduces this.

### 4. Alignment comparison and decision

Development tool, not part of the pipeline: `bin/compare_alignments.py`, run as V2.

- Run nail in full with `-S 1e-2 -E 1000 --ali-out`. Its blocks are simple to parse: `query`, `target`, coordinates, `score`, then the HMM line, match line, target line and posterior line.
- Convert both tools to `hmm_from`, `ali_from` and CIGAR, then compare.
- **Rule (proposed):** nail alone replaces `HMMER_PAIRS` only if all three hold on V2's reads:
  1. GA calls differ on fewer than 0.1% of pairs.
  2. At least 99% of pairs have identical alignments.
  3. Per-family read counts in the final profile match.
- **If it passes:** replace `HMMER_PAIRS` with a script turning nail's `.tbl` and `.ali` into the same domtbl and Parquet, and call GA from nail's score. Change `NAIL_SEED` to a full nail run.
- **If it fails:** keep `HMMER_PAIRS`. Report the differences either way.

The preview above suggests it will fail on the first condition, so build for `HMMER_PAIRS`.

### 5. `SELECTFRAMES`

`modules/local/selectframes/` with `bin/select_frames.py`. Standard library only, using the toolkit container the parser already uses. Logic comes from `measure/orf_filters.py frames` in `kfpfam`.

- **Inputs:** seqkit's ORF FASTA and FGS's `.out`.
- **Outputs:** the kept ORFs with names and sequences unchanged, and a stats JSON (ORFs and residues in and out, reads called, fallback used).
- **Frames:** take FGS's frame at both ends of each call, in forward read coordinates. A frameshifted call therefore marks both frames.
- **Training:** up to 25,000 reads per fold, two folds split by read order, so each ORF is scored by the table it did not train. Pseudocount 1 per 4,096 hexamers.
- **Fallback:** below a minimum number of calls, keep all ORFs and record it. Set the minimum from V3's subsampling, not by guesswork.
- **Check:** a `pytest` file covering the seqkit-to-FGS frame mapping on both strands against real translations, plus the keep rule on a small hand-built case.

### 6. `FRAGGENESCANRS`

`modules/local/fraggenescanrs/`, stock `quay.io/biocontainers/frag_gene_scan_rs:1.1.0--h4349ce8_0`, `-t illumina_5 -w 0`. Publish nothing from it. Use only coordinates from `.out`, never its proteins, because 89% of calls in testing had a repeated codon at a read end.

### 7. Wiring

- **New subworkflow:** `subworkflows/local/nail_hmmer/` takes `orfs` and `pfam_db`, and emits `domtbl`, `alignments` and `versions`.
- **`profile_hmmsearch_pfam/main.nf`:** add the filter branch and engine switch, keeping the current chunking block unchanged under `hmmsearch`. Add an `alignments` emit.
- **Config:** add `conf/modules.config` entries, including resources with retries, and publish `${meta.id}_pfam_alignments.parquet` beside the profile.
- **Parameters and docs:** add them to `nextflow.config` and `nextflow_schema.json`, and update `docs/output.md`, `CITATIONS.md` and `CHANGELOG.md`. Replace the placeholder text in the subworkflow's `meta.yml`.

### 8. nf-test

On the existing mini fixtures and `pfam_mini`:

- Default settings match the existing snapshots, unchanged.
- Filter on: kept ORFs are a named subset of seqkit's.
- `nail`: profile and Parquet are produced, and the Parquet has the expected columns.
- `nail` with the `not_matches` fixture: empty outputs and a clean finish.
- Both options on together.

## Validation

Local runs: one at a time, at most 6 threads. Large runs ship as a Nextflow workflow for the user to run on the HPC; nothing here touches the HPC directly.

| | Test | Pass condition |
|---|---|---|
| **V1** | `hmmer_pairs.py` on every (model, ORF) pair for 10,000 reads, against `hmmsearch -Z 30134 --cut_ga --noali` | Same domains, scores and coordinates; profile identical |
| **V2** | Nail's and HMMER's alignments on 50,000 merged fragments with seqkit ORFs (step 4) | Applies the rule in step 4 |
| **V3** | The same fragments, four arms: current, filter, nail, both. Includes subsampling to set the filter's minimum calls | Compare final `_pfam_coverage.tsv` and stats per family: read counts, breadth, worst families. Your call on acceptable loss |
| **V4** | About 1 M reads through nail and pairs, timed, with peak memory | Cost per read and memory known; decide on chunking |
| **V5** | HPC: other environments and read lengths (soil, marine, eukaryote- or virus-rich, high and low GC; 100, 150, 250 nt), four arms | Does −5 hold across read lengths, since scores sum over ORF length? |

Stop and report if V1 fails, since the design depends on it. V5 decides whether either option could become the default; this branch leaves both off.

## V5 results

Eleven ENA samples through `hpc/pfam-arms` on Codon, Pfam 38.2 in 20 chunks, 2 M merged reads
subsampled per sample. Two runs were needed: in the first, `current` re-drew its subsample while the
other three arms reused a cached draw, so the arms were not comparable. The figures below are from
the second run, where all four arms share one draw, `SEQKIT_SHUFFLE_FASTA` and `SEQKIT_TRANSLATE`
being cached in `filter`, `nail` and `both`. Cost comes from the arms' Nextflow traces, core-hours as
realtime times %cpu, since Codon's trace has no `cpus` column.

| arm | Pfam steps | core-h | wall-h | peak RSS |
|---|---|---|---|---|
| current | `HMMER_HMMSEARCH` | 549.9 | 251.6 | 0.1 GB |
| filter | FGS + `SELECTFRAMES` + `HMMER_HMMSEARCH` | 357.8 | 164.3 | 1.1 GB |
| nail | `NAIL_SEED` + `HMMER_PAIRS` + `MERGEPARQUET` | 26.9 | 8.9 | 16.6 GB |
| both | FGS + `SELECTFRAMES` + `NAIL_SEED` + `HMMER_PAIRS` + `MERGEPARQUET` | 16.7 | 5.3 | 10.0 GB |

`both` is about 33x cheaper than `current` and a third cheaper than `nail` alone: the filter makes the
seeding cheaper too, 13.6 against 22.4 core-h, and takes its peak memory from 16.6 to 10.0 GB. The
filter keeps a median 62% of residues and 63% of ORFs, 96% of reads get a gene call, and the fallback
fired on no sample.

Against `current`, with the draw held fixed, over all eleven samples:

| arm | reads moving, median (range) | families lost | gained | families with 10+ reads below 0.8 |
|---|---|---|---|---|
| filter | 0.18% (0.06-0.33) | 14-54 | 0-2 | 5-24 |
| nail | 2.12% (0.58-2.87) | 19-70 | 0-3 | 12-112 |
| both | 2.26% (0.81-2.91) | 54-105 | 0-3 | 23-118 |

Those are out of 4,600-8,300 families and 2,500-5,400 families with at least 10 reads. So the ORF
filter is close to free, 0.18% of reads for a third off the search; the engine costs about 2% of
reads, one-directional, in line with V2's 1.67% of pair calls; and the two losses are nearly
additive, `both` landing about 0.14 points above `nail`. Both options are cheap enough in accuracy
that cost, not loss, decides them.

**Retraction.** Two earlier readings of this comparison were wrong, and both came from the first
run's mismatched subsamples, which made every arm differ from `current` by 8-12% of reads. That was
neither `PARSEHMMSEARCHCOVERAGE`'s greedy de-overlap flipping near-equal families, as first written,
nor a set of hits that only `current` reported, as written after that. It was a different 2 M reads.
The de-overlap work still stands on its own and says the resolution is rarely contested: in the nail
arms only 2.4% of reads have an overlapping hit of another family at all, the rival is within 2 bits
for 0.4% of reads and within 5 bits for 0.9%, and applying the ORF filter with the engine fixed left
99.7% of 324,247 reads on the same family, 20 of them changing and 997 dropping out. The lesson is
the methodological one: an unseeded subsample produced differences the size of the effects being
measured, in both directions, and looked like signal. Seeded since `8d81c2b`.

**Still open.** V5's own question, whether the filter's -5 threshold holds across read lengths, needs
the per-sample read lengths set against the filter's loss; the spread above, 0.06% to 0.33%, is the
place to look. Separately, `SEQKIT_SHUFFLE_FASTA` peaks at 39 GB, which is `seqkit faidx` indexing
every read in the full FASTA rather than anything to do with the subsample, so it grows with library
size and not with the 2 M cap; `seqkit grep -f` or a one-pass selection would stream instead.

## Order of work

Commit after each step:

1. Step 1, Pfam 38.2.
2. Step 2 and V1, the largest risk.
3. Steps 3 and 4 and V2. This decides nail alone or nail with HMMER.
4. Step 7's engine switch.
5. Steps 5 and 6, then the filter branch of step 7.
6. Step 8.
7. V3 and V4, logged in this file.
8. Ship the V5 workflow, then open one PR against `master` when you ask for it.

## To check while building

- **Biocontainer under Singularity:** `nail` and `mmseqs` run in it under Docker emulation; confirm under Singularity, and record versions in `versions.yml`.
- **`pyhmmer` version:** record which HMMER it bundles, because the pipeline's `hmmsearch` is 3.4. V1 tests the results either way.
- **Read names:** confirm seqkit, nail and FGS handle merged-read names the same way, including any with spaces.

## Depends on

Branch `fix-fetchunzip-checksum` makes `FETCHUNZIP`'s checksum test fire. Until it is merged, step 1's checksum is not enforced.

## Progress log

- 2026-11-06: plan written; test figures above come from `kfpfam`'s `ERR10889147` data.
- 2026-11-06: Pfam set to 38.2, `pyhmmer` set to gathering thresholds, nail-against-HMMER comparison added, at the user's request.
- 2026-11-06: step 1 done. `FETCHUNZIP`'s checksum fix is in PR #29; not yet on this branch.
- 2026-11-06: step 2 done, `HMMER_PAIRS` with nf-test on `pfam_mini`. Run tests as CI does: `NXF_VER=25.10.4 nf-test test ... --profile +docker`.
- 2026-11-06: **V1 passed.** 10,000 reads sampled from merged-fragment seqkit ORFs (62,308 ORFs) × 30,134 models, against `hmmsearch -Z 30134 --cut_ga` (HMMER 3.4, `kfpfam-tools`): all 1,753 domains match in every domtbl column; profile and stats identical; CIGAR, residues, posterior characters and coordinates match `hmmsearch`'s printed alignments, and the CIGAR round trip holds. 67 bytes per domain in Parquet. 284 s wall against 377 s for `hmmsearch` at 6 threads, after fixing a thread contention bug: blocks built across 6 threads ran about 40× slower than a shared block, so they are built in the main thread.
- 2026-11-06: pair step's cost on 3.2 M nail-seed pairs, one thread, about 26.5 core-s: search about 16.7 (MSV about 2, Viterbi about 5, Forward about 6.5, domains about 3), HMM parsing 7, blocks 2, reading input 1.5. Nail seeding is about 400 core-s on the same reads, so compiled changes to HMMER are not worth pursuing unless V4 shows otherwise.
- 2026-11-06: step 3 done, `NAIL_SEED` with nf-test in the biocontainer, under Docker emulation. `-S` is passed through `ext.args` in step 7.
- 2026-11-06: step 4 done. **V2: nail alone fails the rule, so `HMMER_PAIRS` stays.** 50,000 merged fragments with seqkit ORFs (311,680 ORFs), nail 0.7.1 with MMseqs2 18, `-S 1e-2 -E 1000`, against `HMMER_PAIRS` on the same seeds. GA pairs: HMMER 8,448, nail 8,551, both 8,428, nail only 123, HMMER only 20, so calls differ on 1.67% (rule: under 0.1%). Of pairs called by both, 98.39% have identical alignments (rule: at least 99%); on average 99.81% of aligned columns are shared. nail minus HMMER score: median 0.02 bits, 1–99% −0.54 to 0.87. Profile: 103 of about 1,900 families differ, 116 reads in absolute difference out of 8,119 (1.43%), worst PF01473 at 5 against 1. Nail reports one alignment per pair; HMMER has more than one domain on 13 of 8,448 pairs. Time: nail in full 277 s wall, 871 core-s; its seeding stage 189 s wall, then `HMMER_PAIRS` 14 s wall, 30 core-s, 482 MB. Parquet 60 bytes per domain, 508 kB for 50,000 fragments.
- 2026-11-06: step 7's engine switch done: `subworkflows/local/nail_hmmer`, `pfam_search_engine`, `pfam_nail_seed_pvalue`, published `*_pfam_alignments.parquet`. With the default engine, the end-to-end snapshots are unchanged. With `nail`, the end-to-end test on the mini data gives 71 reads in 21 families against 72 in 22. An invalid engine is rejected by the schema. `conf/test.config` had `cpu:` instead of `cpus:` in `resourceLimits`, so CPU limits never applied; fixed. The four older `PROFILE_HMMSEARCH_PFAM` tests and the `FETCHDB` tests fail on `master` too and aren't run by CI: the former have no `ext.args` for seqkit outside the full pipeline's process name, the latter have outdated snapshots. Left alone.
- 2026-11-06: steps 5 and 6 done: `SELECTFRAMES`, `FRAGGENESCANRS`, `pfam_orf_filter` and `pfam_orf_filter_threshold`, with the filter's stats published. Training uses evenly spaced called reads rather than the first ones, because merged reads may not come in random order. `select_frames.py`'s frame logic is checked against real seqkit translations on both strands. `--min-calls` defaults to 5,000 as a placeholder, to be set from V3's subsampling. On the 50,000 merged fragments: 59.5% of residues kept (9.38 M to 5.58 M); 8,436 of 8,448 GA pairs from nail and HMMER on all ORFs are on kept ORFs (99.86%), 3 s and 37 MB for `select_frames.py`. That keeps a larger share than the 54% seen on unmerged reads.
- 2026-11-07: step 8 done. Added a combined filter-and-nail end-to-end test, not yet run because Docker failed (below), and `modules/local/hmmer_pairs/tests/test_hmmer_pairs.py`, which checks the Parquet schema, that CIGAR, residues, posterior characters and coordinates have consistent lengths, and that output doesn't depend on threads (passes).
- 2026-11-07: **V3, 50,000 merged fragments, four arms, final profiles.** Against `current` (`hmmsearch` on all ORFs, 8,195 reads in 1,924 families):

  | Arm | Reads | Families | Read-count difference | Families lost | Of 177 families with 10+ reads, below 0.8 | GA pairs kept |
  |---|---|---|---|---|---|---|
  | filter | 8,186 | 1,924 | 0.11% | 0 | 0 | 99.86% |
  | nail | 8,119 | 1,916 | 0.95% | 8 (15 reads) | 1 (PF00070, 12 to 8) | 99.00% |
  | both | 8,110 | 1,916 | 1.06% | 8 (15 reads) | 1 (PF00070, 12 to 8) | 98.86% |

  No arm found a GA pair that `current` missed, and every shared pair had identical domains, so the arms only lose hits. `hmmsearch` took 2,689 core-s on all ORFs and 1,632 on the filtered ones (61%). With the filter, nail's seeding took 362 core-s and `HMMER_PAIRS` 19, against about 870 and 28 without it in V2.
- 2026-11-07: `--min-calls` set to 1,000 from `kfpfam`'s `measure/filter_min_calls.py`. Training on N called reads and scoring all 50,000: 100 reads keep 66–67% of residues and 99.5–99.6% of GA ORFs; 1,000 keep 57–58% and 99.7–99.8%; 25,000 keep 59% and 99.9%.
- 2026-11-07: **V4, only partly done.** On all 1.04 M reads, `select_frames.py` took 40 s and 410 MB and kept 59.5% of residues. FragGeneScanRs took 41 core-s. Nail's seeding on the 3.9 M filtered ORFs then failed when its MMseqs2 temporary files filled the disk, so nail's temporary files are much larger than its seed file. A second attempt put them inside Docker, whose disk image lives on the same Mac drive. That filled the drive, and Docker's data disk went read-only with I/O errors, so the run stopped. Time, memory and temporary disk at 1–2 M reads now come from V5's traces. Set `NAIL_SEED`'s scratch space accordingly.
- 2026-11-09: `pfam_search_engine = nail` now errors on a chunked Pfam database instead of searching chunk 1. Found while preparing the V5 run on Codon, whose production config serves Pfam as 20 chunks: `NAIL_HMMER` takes `pfam_db.first()`, so the profile would have been about 1/20 complete with nothing in the log to say so, and `-Z num_models` wrong per chunk. The guard is in `PROFILE_HMMSEARCH_PFAM`, where the engine is chosen, with an nf-test (`tests/chunked_pfam.config`). Chunked nail would need per-chunk seeding, a domtbl concatenation and a Parquet merge; not worth it while `hmmsearch` is the default engine. All six end-to-end tests pass, including the filter-and-nail test that Docker had blocked.
- 2026-11-07: V5 shipped in `hpc/pfam-arms`: `run.sh` (sbatch), `site.config` (Slurm, Singularity), nine ENA samples in `samplesheet.csv`, `bin/compare_profiles.py` and `cost.py`. Not run.
- 2026-11-10: chunked Pfam supported by the `nail` engine, replacing yesterday's refusal, so the V5 run can use Codon's 20-chunk Pfam. `NAIL_HMMER` fans each sample out over the sorted chunks (`hmm_chunk` in the meta keeps the seeds joined to their own chunk), then `COMBINEHMMSEARCHTBL` concatenates the domtbls and a new `MERGEPARQUET` merges the alignments, sorted and written as `HMMER_PAIRS` does. Both run for a single chunk too, so there is one code path; the Parquet is published from `MERGEPARQUET` instead of `HMMER_PAIRS`. `-Z` is still `databases.pfam.variables.num_models`, the whole database either way. `NAIL_SEED` and `HMMER_PAIRS` now take the HMM in their input tuple, since it varies per task.
- 2026-11-10: verification of chunked nail on this Mac was partial, because every `nextflow` run here now dies with `Unable to create serializer ... ThreadPoolExecutor` (kryo against `nextflow.Session`) part-way through, on `master` as well as on this branch, under Nextflow 26.04.1 and JDK 23. What was confirmed on the two-chunk test config: `NAIL_SEED` and `HMMER_PAIRS` each run once per chunk, and the concatenated domtbl is 152 lines from two 76-line chunks. `MERGEPARQUET` is killed by that abort mid-run, so its script was run by hand in its container on the two chunk Parquets: 76 + 76 = 152 rows, schema preserved, sorted by `pfam`. The nf-test `workflow.trace.succeeded()` assertions use `it.name`, since `WorkflowTask` has no `process` property. Re-run the suite where Nextflow works before trusting it green.
- 2026-11-12: B11's results read. Cost over 11 samples: `current` 540 cpu-h in `HMMER_HMMSEARCH`, `filter` 348, `nail` 26 (21 seeding, 4.5 pairs, wall 8 h against 245 h), so nail is about 21x cheaper and the filter alone saves 36%. Peak RSS: `NAIL_SEED` 16.6 GB, `HMMER_PAIRS` 7.7 GB, and `SEQKIT_SHUFFLE_FASTA` 43.6 GB, the largest step in the run. Chunked nail behaved: 220 `NAIL_SEED` and 220 `HMMER_PAIRS` tasks (11 samples x 20 chunks) and one merged Parquet per sample, 331 MB in all. The profile comparison is confounded, though: `SEQKIT_SHUFFLE_FASTA` subsampled with unseeded `shuf`, and 10 of 11 samples re-ran it in the nail arm, so those profiles are a different draw. On SRR32629731, the one sample whose subsample and translation were cached, nail and the filter differ from `current` by 11.1% and 11.0% of reads but from each other by only 0.93%. Totals barely move (326,759 / 325,526 / 324,247 reads mapped) and `read_hit_count` is about equal to `reads_mapped` in every arm, so `PARSEHMMSEARCHCOVERAGE`'s greedy de-overlap keeps roughly one family per read, and pruning the reported hit set at all — by dropping ORFs or by scoring only nail's candidates — changes the winner for about a ninth of reads. That sensitivity is in the parser, not in either arm. The `both` arm failed and is being re-run.
- 2026-11-12: subsample seeded, so the arms can be compared. Verified by running the module twice through Docker on 1,000 reads: seed 23 gives the same 50 reads both times (md5 500d163b), `ext.seed = 7` a different 50 (fbccea5f). Re-running the arms after this lands will re-execute everything from the subsample down.
- 2026-11-12: the `both` arm finished and V5's results are written up above. Headline: `both` is about 31x cheaper than `current` for 0.1-1.8 points more read movement than `filter` alone, the engine's own effect is about 2%, and most of the 11% against `current` traces to the de-overlap in `PARSEHMMSEARCHCOVERAGE` rather than to either feature.
- 2026-11-13: the arms re-run on one shared subsample, and V5's section rewritten around it. The filter costs 0.18% of reads at the median, the engine 2.12%, both together 2.26%, against `current` over eleven samples; cost is 550, 358, 27 and 17 core-hours. The 8-12% seen in the first run was the mismatched draws, not the parser and not a `current`-only hit set; both earlier readings are retracted in place. Chunking adjusted in `hpc/pfam-arms` from the per-task figures at the same time: `hmmsearch` tasks ran 67 minutes at the median on 2.25 of 4 cores with the reads never split, so `hmmsearch_chunksize` drops to 256 MB and `cpus` to 2, while nail's 20 tasks of 1.6 minutes per sample each rebuilt their MMseqs2 query database, so the nail arms now search an unchunked Pfam built by `run.sh` from the site's chunks, one `NAIL_SEED` and one `HMMER_PAIRS` per sample.
- 2026-11-19: **kfpfam added as a third engine, and the benchmark grew to seven arms**, at the user's request: the pipeline had only ever compared `hmmsearch` with nail, while kfpfam's own README has the three-way numbers locally (1.07 M reads, its 50,000-read truth): Mode A at 1.20 core-ms per read and .999 recall at GA, Mode B at 0.37 and .987, nail's defaults at 1.24 and .987. `pfam_search_engine = kfpfam` runs `filter join` once per sample over the whole index and converts its (ORF, model index, diagonal, gapped) rows into nail's seed format through `prep/models.tsv`'s NAME column, so `HMMER_PAIRS` is unchanged; `subworkflows/local/nail_hmmer` became `seeded_hmmer`, which both engines share. The image, prep directory and one index per mode are built once by `hpc/pfam-arms/kfpfam_setup.sh` (kfpfam is published nowhere and had never run on Codon), with `docker/filter.Dockerfile` added to kfpfam for that. New arms: `nail_loose` at `-S 0.1`, `kfpfam_a` and `kfpfam_b`. Checked by previewing both engines, and the index-to-name conversion by hand on a synthetic `models.tsv`.
- 2026-11-19: `nar-review-response/hpc/B11_pfam_arms/arms.config`'s unchunked-Pfam block **had never taken effect**: assigning `params.databases.pfam` merges into `codon.config`'s map, so `chunked: true` and the twenty chunk entries stayed behind. That is why the first B11 run shows 220 `HMMER_HMMSEARCH` and 220 `NAIL_SEED` tasks, 20 per sample, and why the nail arms first met the chunked-database guard. Fixed with `chunked: false`, so the arms now search one whole Pfam and nail and kfpfam seed once per sample. The `--tmp=300G` request that blocked submission on 2026-11-10 is gone, and the task resources there come from B11's measurements. This repo's `hpc/pfam-arms` samplesheet was a stale nine-run DRR set; it now matches the live benchmark, six Zymo standards plus five biome runs, with the two 250 nt runs kept for the read-length question.
- 2026-11-19: kfpfam image built, `kfpfam-filter:latest`, amd64, 60 MB, from kfpfam at `cfefbe4`; a `.dockerignore` was needed first, since that checkout's untracked `data/` is 30 GB of build context. Exported with `docker save` to `nar-review-response/hpc/B11_pfam_arms/kfpfam/kfpfam-filter.tar` (57 MB), because neither Singularity nor Apptainer is installed here, and `kfpfam_setup.sh` now converts that archive with `docker-archive://` on Codon, which has no Docker either. A smoke test through the image on the mini Pfam fixture (34 models, 853 ORFs) found a real bug: the filter's reader asserts one line per sequence while `seqkit translate` wraps at 60 columns, so every real sample would have panicked; `KFPFAM_FILTER` now flattens the ORFs before calling it and leaves the wrapped file for `HMMER_PAIRS`. End to end after the fix: 197 candidate pairs, converted to seeds naming models by their NAME field, and `hmmer_pairs.py` on those gave a 65-row domtbl and a Parquet.
- 2026-11-19: alignment publishing checked, at the user's request, with a follow-up in mind that anchors reads to MGnify proteins through their Pfam alignments. `MERGEPARQUET` publishes `*_pfam_alignments.parquet` per sample for both seeded engines, zstd level 9, dictionary-encoded families, sorted by family then read; B11 measured 14-52 MB per sample at the 2 M cap, 331 MB for eleven samples, about 60 bytes a domain. Added `read_from` and `read_to` so the anchor needs no frame arithmetic. Writing them exposed a defect in `hmmer_domtbl_parse_coverage.py`: its nucleotide mapping adds the frame offset a second time on top of the ORF's `begin`, and places a reverse frame one nucleotide low, so 13 of 65 fixture domains fell outside their own ORF. Harmless where it is used, since the shifts cancel when it compares hits on one read, but it means the parser's coordinates are not read coordinates; `read_span` computes them properly and `test_hmmer_pairs.py` checks both strands, all three phases, and that a whole ORF maps back to `begin..end`. The `hmmsearch` arms still publish no alignments, since the search runs `--noali` and only its domain table survives: that would need a domtbl-to-Parquet step, with coordinates but no CIGAR.
- 2026-11-19: Pfam chunking made per arm, at the user's request, since the engines want opposite things: `hmmsearch` parallelises over models, so the `current` and `filter` arms keep Codon's 20-way chunked database, while nail and kfpfam read one whole HMM file and now get an unchunked one. It was wrong both ways beforehand. In the live benchmark the unchunked block sat in `arms.config`, which every arm uses, so my `chunked: false` fix had just moved the mistake: all seven arms would have searched one whole Pfam. It now lives in `nar-review-response/hpc/B11_pfam_arms/unchunked_pfam.config`, which `launch.sh` passes to the seeded arms only, and the hmmsearch arms are back on the database they were first run with, so their cached tasks resume. In this repo's `hpc/pfam-arms/run.sh` the two kfpfam arms had been left without `-c unchunked_pfam.config` at all, so they would have seeded once per chunk; they have it now.
- 2026-11-19: subsampling turned off for the benchmark, at the user's request (`--hmmsearch_subsampling -1` in both launchers), and the arms reordered so the seeded ones run first: whichever arm runs first pays for the shared download, QC, decontamination and merge, and the cheap arms' real numbers should decide whether the hmmsearch arms are worth running at full depth. The eleven libraries hold about 440 M read pairs against the 22 M B11 searched, so on B11's per-read cost `current` is of the order of 10,000 core-hours and `filter` 7,000, against a few hundred for each seeded arm. Resources raised accordingly, and two of them are now unmeasured rather than measured: `NAIL_SEED` is one task per sample over a whole library (120 GB, 48 h, and its MMseqs2 temporaries are what filled a disk in V4), and `HMMER_PAIRS` holds every reported domain in memory before writing the Parquet, as `MERGEPARQUET` does to sort it. If either is killed for memory the fix is to stream the Parquet, writing row groups as the search goes and sorting on read, not to raise the request again. `SEQKIT_SHUFFLE_FASTA` drops out of the run entirely, which also removes the 39 GB step.
- 2026-11-19: the kfpfam arms failed on Codon for two reasons, both fixed. Nextflow's task wrapper runs `ps` inside the container (`nxf_trace_linux`, `command-trace.txt`) and **exits 1** when it is absent, so every `KFPFAM_FILTER` task died on `debian:trixie-slim`'s lack of `procps`; the image now installs it, and `kfpfam_setup.sh` rebuilds the `.sif` when the exported tar is newer so a stale one cannot survive. Separately, `--pfam_search_engine kfpfam` was rejected by the schema because the shared asset clone at `/homes/emgpr/.nextflow/assets` was behind the branch: `-r <branch>` and `-resume` do not refetch, so `nextflow pull` is needed before a launch that depends on new code. Worth knowing for local testing: Nextflow 26.10 rejects every integer passed on the command line as a string, including `--hmmsearch_subsampling=-1`, while Codon's version accepts it — put integers in a config file to test here.
- 2026-11-19: **read chunking for the seeded engines**, at the user's request, 4 M reads per chunk (`pfam_seeded_chunk_reads`, 0 to disable). The seeded engines read one whole HMM file, so they have no database chunking to spread a sample over, and with subsampling off that made translation, seeding and HMMER one task each per 40 M-read library. The split is on the reads, before translation, not on the ORFs: 4 M is then exact rather than a guess at ORFs per read, and it also chunks `SEQKIT_TRANSLATE`, which `conf/modules.config` allows one hour and 1 GB — sized for the 2 M cap, and a certain failure at full depth for every arm, the hmmsearch ones included, which still translate a whole library in one task. Regrouping closes over both axes, `hmm_chunks * read_chunks`, and the per-chunk outputs carry zero-padded `r%04d_h%04d` prefixes so they neither collide when staged side by side nor concatenate in finish order. Verified on a standalone harness of the real operators (7 reads at 2 per chunk x 3 Pfam chunks = 12 tasks, 1 read = 3, groups closing at exactly those sizes), because a full `-stub-run` cannot complete: `rawr.nf:159` parses FASTP's stub JSON, which is empty, on `master` too. Expected cost per sample for a kfpfam arm, extrapolated from B11's nail numbers and so not measured: ten chunks, about 50 core-hours, 1-2 h wall once its chunks are running. `MERGEPARQUET` is the one step chunking makes heavier — ten inputs per sample instead of one, the same rows, still sorted in memory — so it is now the likeliest to need the streaming fix.
- 2026-11-19: care with the live `nar-review-response/hpc/B11_pfam_arms/arms.config`: it is **not** interchangeable with this repo's `hpc/pfam-arms/site.config`, and copying one over the other loses things silently. `arms.config` layers on `codon.config`, so it carries resource overrides only, plus `cleanup = false` — without which the first arm to finish deletes the work directory every other arm resumes from, since this repo's `nextflow.config` sets `cleanup = true` — a `publishDir` for `BBMAP_REFORMAT`'s merged reads, fully-qualified process selectors, and no `hmmsearch_chunksize`, which it defers to `codon.config`. `site.config` is standalone and declares its own `executor`, `cleanup`, error strategy and chunk size. Keep them separate.
- 2026-11-23: **V6, the first full-depth run: kfpfam A and B against nail on all eleven runs**, written up in `nar-review-response/response/B11_v6_kfpfam_arms.md` from `scripts/b11v6_kfpfam_arms.py`. The three seeded arms are complete and searched identical reads; `current` is not (1,679 of about 2,910 `HMMER_HMMSEARCH` tasks, no `PARSEHMMSEARCHCOVERAGE` at all), so the profiles under `results_main_current` are still the 2 M-subsample ones and are usable only as a depth reference. Uncapping the search is worth more than the choice of prefilter: recall over all truth families goes from 0.735 at 2 M to 0.93. kfpfam tracks nail closely (Bray-Curtis 0.013, Spearman 0.999 over about 10,000 shared families, Mode A annotating 0.9% more reads and Mode B 1.3% fewer) and **recovers exactly the families nail was losing** -- PF01391 Collagen at 8.0x nail's reads, with PF05594, PF01521 and the rest of nail's worst list at 1.4-1.8x -- which answers the question the previous report left open. Against that, **183 families that nail finds are identically zero in both kfpfam modes**, 0.78% of reads: Mode A is strictly more sensitive than Mode B, so these models are absent from the k-mer index rather than missed by it, and dropping them from truth closes the whole accuracy gap (Mode A then beats nail, +0.002 recall, +0.003 rho). Check them against `prep/models.tsv` on Codon before anything else. Cost upends kfpfam's README: Mode A took 2,153 cpu-h against nail's 287 and Mode B 585, so neither mode is cheaper than nail here, though their HMMER stages behave as designed (17 cpu-h for Mode B against nail's 81). Everything outside Pfam profiling is 22 cpu-h in total, 0.15% of `current`'s projected 14,600. Memory is the other way round: `NAIL_SEED` needs 68 GB at the median and 122 GB at the worst per 4 M-read chunk and **98 of 229 attempts were OOM-killed at 120 GB**, while `KFPFAM_FILTER` holds 5 GB flat. So the resource cut of 2026-11-19 was wrong -- `NAIL_SEED` is back to 120 GB in both configs, `KFPFAM_FILTER` down to 16 GB, and both now carry the measured figures rather than estimates.
- 2026-11-23: **`speed_test.sh` on Codon, 800,000 ORFs of SRR23926920, Mode A, 8 threads:** kfpfam `f40154c`'s AVX2 kernel 347 core-s, SSE4.1 374, the old path every B11 arm ran on 1,053, identical pairs. Mode A's 2,153 cpu-h would be about 710, still 2.5x nail's 287, with the gapped stage (177 core-s) now the largest. kfpfam then got a SIMD bound ahead of that stage and a 32-lane AVX2 kernel, both lossless (`kfpfam/docs/scaling.md`); not timed on Codon yet. The 32 lanes want larger blocks, so `kfpfam_block` (400,000) is to be set from `speed_test.sh`'s block sweep.
- 2026-11-23: **the kfpfam arms searched the wrong Pfam, which explains the missing families and part of the cost.** `diagnose_183.sh` on Codon: `kfpfam/prep` was built at 14:56 on 2026-11-19 from Pfam 39.0 (`current_release`, 53,275 models; the file is the same 658,389,433 bytes), before the arms had downloaded 38.2 (15:07), while `HMMER_PAIRS` read 38.2's 30,134. Seeds join to models by name and an unknown name was skipped in silence, so renamed families got nothing (PF00577 is B-barrel_usher there and Usher in 38.2). Of the 183, 20 are not in that prep at all, 144 are under another name and 19 under the same name; those 19 are not explained by the join, and need the re-run to check. 796,262 of 2,226,587 seed rows (36%) named models 38.2 lacks. The search cost about 1.7 times what 38.2 would (3,420 hits per ORF residue; 1.77 times the models), so V6's kfpfam costs are too high and its accuracy is that of a mismatched index: **V6's kfpfam columns should not be used**, and both arms need re-running after `kfpfam_setup.sh` rebuilds. Guards added: `kfpfam_setup.sh` refuses a file with the wrong number of models and rebuilds a prep that was not made from `PFAM_HMM` (by checksum), `KFPFAM_FILTER` stops if prep and `num_models` disagree, and `hmmer_pairs.py` raises on a seeded name the whole database lacks. The faster kfpfam build timed since is 157 core-s against `f40154c`'s 361 on the same node and 800,000 ORFs, and 2.3x on this job's 400,000 too (75 against about 173), but `ac87746` writes one pair fewer, which `kfpfam/docs/scaling.md` explains and has a fix for. Projected, not measured: that build against 38.2 would put Mode A near 180 cpu-h over the eleven runs, against nail's 287 (2,153 x 75 / 526, divided by 1.7).
- 2026-11-23: **moved to Pfam 39.0** (53,275 models), the user's decision, after the mix-up above. `nextflow.config` pins `releases/Pfam39.0` by md5 with `num_models` 53275. The pipeline's download cache is found by name, not release, so it must be refreshed with `--force_download_dbs`; in B11 every arm reads `$B11_DIR/pfam/39.0` by `local_path` instead (`pfam_setup.sh`, whole for the seeded arms, 20 chunks for the hmmsearch ones), and results go under `pfam39.0/`, apart from the 38.2 runs. The `hmmer_pairs.py` check is now `--whole`, passed when Pfam is one file: by model count it would have let exactly this case, a 39.0 prep with a cached 38.2, through as a chunk. Consequences: nail and the hmmsearch arms search again too, production's `chunked-20` has to move with the `-Z` (its release is unchecked), and B9's truth is 38.0 until `truth_refs_job.sh` is re-run. The speed test's index on Codon was 39.0 all along, so its figures stand for this database: 157 core-s against the 1,053 of the build V6 ran puts Mode A at about 320 cpu-h over the eleven runs (projected). nail's 287 was against 38.2, so the two are not comparable until nail has run on 39.0.
- 2026-11-23: **ready for the 39.0 benchmark of `kfpfam_a`, `kfpfam_b` and `nail`.** kfpfam `7b3bf57` in the last `speed_test.sh` (EPYC 9555, 8 threads, 800,000 ORFs, 39.0): Mode A 121 core-s against 354 for `f40154c` and 1,053 for the kernel V6 ran, every kernel's pairs identical; Mode B 52-56. Projected from V6's 2,153 and 585 cpu-h: about 250 and 110. HMMER on kfpfam's candidates will cost more than V6's 64 and 17, as that run dropped 36% of its seed rows. nail's 287 was against 38.2 and may grow with the models (1.77x), so `NAIL_SEED` asks 240 GB and `HMMER_PAIRS` 96 in B11. `cost.py` could not read Codon's traces (no `process`, `tag` or `cpus` column) and left `KFPFAM_FILTER` and `MERGEPARQUET` out; it reproduces V6's 286.8, 2,152.9 and 585.0 cpu-h now. B11's `launch.sh` takes any arms and chains them, checks everything before submitting, and pulls the branch; `collect.sh` no longer needs the hmmsearch arms. `profile_hmmsearch_pfam`'s nail tests pass (Nextflow 26.04 needs `NXF_SYNTAX_PARSER=v1` for nf-test); its four `hmmsearch` snapshots fail on `versions.yml` hashes, with and without these changes.
- 2026-11-25: **`pfam_weak_evalue`, off by default**, at the user's request, for the MGnify annotation index (kmer-functional-profiler, `workflows/pfam-anchor/lookup/ANNOTATION.md`): only 41% of 50-residue windows of KO proteins align at the gathering thresholds. The seeded engines' Parquet gains domains under them, where the ORF's and the domain's E-values are at most the parameter, and the columns `above_ga` and `pfam_version`. E is at the pipeline's Z, so 0.025 on 37.4 is the P-value of 1e-6 that was tested. Each candidate is searched without a cut first, and at the gathering thresholds only the ORFs that search reports. On 12,000 windows with kfpfam's candidates (k4 T7.5) and Pfam 38.2, the domtbl is the same byte for byte, the 4,115 GA rows too, 3,556 weak rows are added, all of them `hmmsearch`'s, and HMMER takes 7 core-s against 6, model loading included. KO recall of all such windows goes from .360 to .468, at precision .957 and .949. Tested by `test_hmmer_pairs.py` and by running `SEEDED_HMMER` with Nextflow, natively and offline; **not** by nf-test, whose plugin is not on this laptop, nor with nail or in a container. The new process input changes the task hashes, so a resumed run searches again.
- 2026-11-25: **`pfam_lookup`, off by default**, at the user's request: `PFAM_LOOKUP` profiles each sample's alignments against an index of MGnify's proteins (kmer-functional-profiler's `workflows/lookup-index`, its `bundle.py profile`), and publishes `*_pfam_lookup.*`. Off-by-default, and an error with the `hmmsearch` engine, without `pfam_alignments`, or without the bundle and image. `hpc/pfam-arms/lookup_setup.sh` makes the image's `.sif` and links the bundle; the kfpfam arms use them when both are there. Checked: the module's real command on a fixture bundle (`modules/local/pfam_lookup/tests/test_pfam_lookup.sh`, 18 of 18 reads called), its nf-test as a stub, the errors in `-preview`, and a stub run from `SEEDED_HMMER` to the published files under Nextflow 25.10.4. Not checked: inside the image (not built here, for disk), and a whole pipeline run (`FETCHUNZIP`'s stub fails on `rfam`, before this). The bundle is of Pfam 39.0, as the arms (the user's decision the same day; 38.2 was planned). Arms last run before `pfam_weak_evalue` search again, as `KFPFAM_HMMER`'s command changed; `kmer-functional-profiler/hpc/kfp-ablations/lookup-profile` profiles alignments that exist already. Since then also checked inside the image (amd64, emulated): `PFAM_LOOKUP` calls the fixture's 18 reads.
- 2026-11-27: **kfpfam top with the MGnify lookup is the standard run**, at the user's request. The pipeline's kfpfam defaults are top's (bar 15, gate 14.5, HMMER at 21.25), `KFPFAM_HMMER` exports `FILTER_GAPPED_FILE` when the prep has `gapped.txt` and runs with `FILTER_BAND=15 FILTER_DIAGONALS=any` (`conf/modules.config`; `weefams.config` no longer sets them). `hpc/pfam-arms/kfpfam_setup.sh` builds top (theta 9, seed 10110001, gapped bars at slack 0, hm-annotate's PREP) from the lookup bundle's own Pfam (38.1) into `kfpfam/`, keeping that `Pfam-A.hmm`; `run.sh`'s `kfpfam_top` arm searches it with `--pfam_lookup`, replacing `kfpfam_a` and `kfpfam_b`; `RUN_ARMS=kfpfam_top` runs it alone, and the comparison is skipped without `current`. Checked: top's build steps and `filter join` with those settings on `tests/fixtures/pfam_lookup` (kfpfam 9833cf6's binary, natively: 18 of 18 ORFs paired; a missing gapped file is an error, so it is read), `nextflow config` with and without `weefams.config`, `bash -n`. Not checked: inside the image, on Codon, or a pipeline run.
