# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

`cnr-ibba/nf-blobpurge` is an nf-core-style Nextflow pipeline that runs **after** [sanger-tol/blobtoolkit](https://github.com/sanger-tol/blobtoolkit) has already classified a short-read de novo genome assembly. Given an existing assembly + BlobDir pair per sample, it:

1. **`BTK_FILTER`** — removes contaminant contigs via `blobtools filter` against the existing BlobDir, with a programmatic span-conservation check (`bin/check_span.py`) that fails the process if `span(filtered) + span(excluded) != span(original)` within `--span_tolerance`.
2. **`READ_COVERAGE`** — produces a purge_dups coverage stat/base-cov pair: subsets an existing reads CRAM (already aligned upstream), or maps `reads_r1`/`reads_r2` fresh with `bwa-mem2` onto the *filtered* assembly.
3. **`PURGE_DUPS`** — always runs: `split_fa` → self `minimap2` → `calcuts` (+ `hist_plot.py`) → `purge_dups` → `get_seqs`.
4. **`PURGE_HAPLOTIGS`** — optional (`--run_purge_haplotigs`, default `true`), independent cross-check on the *same* coverage input, not a second pass over purge_dups' output. Cutoffs are auto-estimated from the sample's own bimodal depth distribution (`bin/estimate_purgehaplotigs_cutoffs.py`) instead of chosen by eye.
5. **`BUSCO_COMPARE`** — runs BUSCO (never `--auto-lineage`) on every assembly stage (`raw`, `filtered`, `purge_dups`, `purge_haplotigs`) × every requested lineage, then builds one comparative table per sample (`bin/compare_busco.py`).
6. **`BLOBPURGE_REPORT`** — one per-sample HTML report (`bin/generate_report.py`) tying span, BUSCO duplication, the optional GenomeScope2 comparison, and the purge_dups/purge_haplotigs size-disagreement check together into an explicit verdict. Every number in the report is read from an upstream process's own output JSON — nothing is hand-entered, and missing optional inputs are reported as explicitly "not available."
7. **MultiQC** — aggregate QC/versions report across samples.

Nothing about organism/taxon/lineage is hardcoded: `--exclude_taxa` and `--busco_lineages` are **mandatory** parameters with no default; the pipeline fails immediately and explicitly if either is missing.

Known caveats the pipeline surfaces on purpose (both to stderr and in the report, via a `ch_caveats` channel threaded through to `BLOBPURGE_REPORT`), rather than papering over:
- purge_dups coverage uses `ngscstat` (the Illumina/short-read path), which is less exercised upstream than the PacBio `pbcstat` path.
- `purge_haplotigs` is designed around long-read coverage; here it typically runs on the same short-read coverage as purge_dups.

## Commands

### Run the pipeline
```bash
nextflow run . -profile test,docker --outdir <OUTDIR>
```
The `test` profile uses the synthetic fixture under `assets/test/` (see `assets/test/generate_test_data.py`) — a 6-contig toy assembly with one intentional near-duplicate contig pair and two "contaminant" contigs, plus matching reads and a hand-built BlobDir. It is a plumbing fixture, not biological data. `conf/test.config` regenerates an absolute-path samplesheet at config-load time (nf-schema resolves samplesheet paths against the launch dir, not the CSV's location).

Real runs require the mandatory params:
```bash
nextflow run . -profile docker \
   --input samplesheet.csv --outdir <OUTDIR> \
   --exclude_taxa "Pseudomonadota,Bacteroidota,Actinomycetota" \
   --busco_lineages "chlorophyta_odb12,viridiplantae_odb12"
```

### Test
```bash
nf-test test --tag test --profile +docker --verbose
```
Update snapshots after intentional output changes:
```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```
There is a single top-level pipeline test, `tests/default.nf.test`, which runs `-profile test` end-to-end and snapshots stable output paths/contents plus the collated versions file. Module/subworkflow-level `.nf.test` files are not currently used in this repo — only nf-core/modules vendored tests exist (ignored by `nf-test.config`).

### Lint
```bash
nf-core pipelines lint .
```
Config obeys `.nf-core.yml` (`repository_type: pipeline`, `is_nfcore: false`, `org: cnr-ibba`). Pre-commit config (`.pre-commit-config.yaml`) and Prettier (`.prettierrc.yml`) are also present.

### Update the parameter schema
After changing `nextflow.config` params, run:
```bash
nf-core pipelines schema build
```

## Architecture

`main.nf` wraps three phases: `PIPELINE_INITIALISATION` (validates/reads the samplesheet) → `CNRIBBA_BLOBPURGE` (calls the `BLOBPURGE` workflow in `workflows/nf-blobpurge.nf`) → `PIPELINE_COMPLETION` (email/notifications). This split is the standard nf-core template shape; the actual pipeline logic lives entirely in `workflows/nf-blobpurge.nf` and the subworkflows it calls.

### Channel shape driving everything
The samplesheet channel is `[ meta, assembly, blobdir, reads_cram, reads_r1, reads_r2 ]`. `meta.has_cram` (set during sample-sheet parsing) drives the CRAM-vs-FASTQ branch inside `READ_COVERAGE`.

A `ch_stage_fastas` channel of `[ meta, stage, fasta ]` accumulates one entry per pipeline stage (`raw`, `filtered`, `purge_dups`, and — if enabled — `purge_haplotigs`) and is the single input both `ASSEMBLY_STATS` and `BUSCO_COMPARE` fan out over. Adding a new assembly-producing stage means mixing a new `[ meta, '<stage_name>', fasta ]` entry into this channel — stats and BUSCO comparison follow automatically.

A `ch_caveats` channel of `[ meta, caveat_string ]` is mixed from any subworkflow that wants to flag a methodological caveat (currently `READ_COVERAGE`/`PURGEDUPS_NGSCSTAT` and `PURGE_HAPLOTIGS`); it's grouped per sample and joined into the final report inputs, so caveats always reach the HTML report rather than being silently logged only to stderr.

### Subworkflow layout (`subworkflows/local/`)
- `read_coverage.nf` — CRAM-subset vs. fresh bwa-mem2 branch, converges on a common `ngscstat` step (needs a name-sorted BAM; coordinate-sorted/indexed BAM is kept separately for `purge_haplotigs`).
- `purge_dups.nf` — linear purge_dups toolchain.
- `purge_haplotigs.nf` — parallel/independent cross-check, reuses the BAM from `read_coverage`.
- `busco_compare.nf` — combines `ch_stage_fastas` with every `--busco_lineages` entry (`Channel.of(*params.busco_lineages.tokenize(','))`), runs BUSCO per (stage, lineage), then one `COMPARE_BUSCO` per sample across all its summaries.
- `utils_nfcore_nf-blobpurge_pipeline/` — nf-core template boilerplate (samplesheet validation, citations text, completion emails). Don't hand-edit generated nf-core template subworkflows under `subworkflows/nf-core/` — those are vendored from nf-core/modules via `modules.json` and updated with `nf-core pipelines modules update` / `nf-core pipelines sync`, not by direct editing.

### Modules (`modules/local/`)
All pipeline-specific logic is a local module (no nf-core remote modules besides `modules/nf-core/multiqc`). Each module directory has `main.nf` (+ `environment.yml` for conda). Notable ones:
- `btk_filter` — wraps `blobtools filter`; auto-detects `--taxrule` from the BlobDir's `meta.json` via `bin/detect_taxrule.py` when not given explicitly, and fails with a pointed error (naming the exact field/taxrule mismatch to check) rather than passing through `blobtools`' raw failure.
- `blobpurge_report` — invokes `bin/generate_report.py`, the final verdict-producing step.
- `purgedups_*` / `purgehaplotigs_*` — one module per toolchain step (mirrors the subworkflow's linear structure), not one big monolithic module.

### `bin/` scripts
Standalone Python scripts invoked from module `script:` blocks, not imported as a library — treat each as the source of truth for the file format it reads/writes (e.g. `assembly_stats.py` emits `<sample>.<stage>.stats.json`, consumed downstream by `generate_report.py` via the naming convention `STATS_FILENAME_RE = ^.+\.(raw|filtered|purge_dups|purge_haplotigs)\.stats\.json$` — keep the two in sync if you rename a stage).

### Versions collection
Uses Nextflow's `channel.topic("versions")` pattern (newer nf-core style) in addition to the traditional `ch_versions` channel mixed through every subworkflow — see the topic-branching/grouping logic at the end of `workflows/nf-blobpurge.nf` before assuming only `ch_versions` needs updating when adding a process.

## Conventions to follow when extending the pipeline

- `--exclude_taxa` and `--busco_lineages` must stay mandatory with no default — this is a deliberate design choice (results must never silently apply an organism-specific default).
- New stage-producing steps: mix `[ meta, '<stage>', fasta ]` into `ch_stage_fastas` in `workflows/nf-blobpurge.nf`; nothing else needs to change for stats/BUSCO to pick it up.
- New caveat-worthy limitations: mix `[ meta, caveat_string ]` into `ch_caveats` from the producing subworkflow rather than only logging to stderr.
- New params go in `nextflow.config` under `params {}` with an explicit default (or `null` if mandatory), then `nf-core pipelines schema build` to sync `nextflow_schema.json`.
- Follow the nf-core channel-naming convention already in use: `ch_<previousprocess>_for_<nextprocess>` for intermediate/terminal channels.
- Resource defaults belong in `conf/base.config` under `withLabel:` selectors (standard nf-core process labels: `process_low`/`process_medium`/`process_high`, etc.), not hardcoded per-module.
- Update `docs/usage.md` / `docs/output.md` and `CITATIONS.md` alongside any new mandatory param, output file, or tool dependency — these are hand-maintained, not generated.
