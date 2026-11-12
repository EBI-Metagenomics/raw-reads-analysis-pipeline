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
