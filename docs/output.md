# cnr-ibba/nf-blobpurge: Output

## Introduction

This document describes the output produced by the pipeline for each sample in the samplesheet. All paths below are relative to the top-level `--outdir` results directory.

## Pipeline overview

The pipeline processes each sample through the following stages:

- [BTK_FILTER](#btk_filter) - contamination removal from an existing BlobDir, with a programmatic span-conservation check
- [Read coverage](#read-coverage) - CRAM subsetting or bwa-mem2 mapping, feeding purge_dups' `ngscstat`
- [purge_dups](#purge_dups) - haplotig purging
- [purge_haplotigs](#purge_haplotigs) - independent cross-check of purge_dups (optional, on by default)
- [Assembly stats](#assembly-stats) - span/N50/contig-count at every stage, machine-generated
- [BUSCO comparison](#busco-comparison) - comparative completeness/duplication across stages and lineages
- [blobpurge report](#blobpurge-report) - the final per-sample verdict
- [MultiQC](#multiqc) - aggregate report across samples
- [Pipeline information](#pipeline-information) - execution metrics and software versions

### BTK_FILTER

<details markdown="1">
<summary>Output files</summary>

- `btk/`
  - `<sample_id>.filtered.fasta`: assembly with `--exclude_taxa` contigs removed.
  - `retained_ids.txt` / `excluded_ids.txt`: contig IDs kept/removed.
  - `<sample_id>.span_check.json`: the programmatic check that span(filtered) + span(excluded, from the BlobDir) == span(original) within `--span_tolerance`. The process fails if this does not hold.
  - `<sample_id>.span_check_mqc.json`: the same span check, reformatted as a MultiQC custom-content table.
  - `<sample_id>.filter_summary.json`: raw `blobtools filter --summary STDOUT` output.

</details>

Runs `blobtools filter` against the existing BlobDir named in the samplesheet, using `--taxon_field`/`--exclude_taxa` (and `--taxrule`, deduced from the BlobDir's `meta.json` if not given).

### Read coverage

<details markdown="1">
<summary>Output files</summary>

- `samtools/` or `bwamem2/`: sorted/indexed BAM used for coverage (CRAM subset, or fresh bwa-mem2 mapping onto the filtered assembly).
- `purgedups/<sample_id>.ngscstat.stat`, `<sample_id>.ngscstat.base.cov`: purge_dups Illumina coverage stats (`ngscstat`).

</details>

### purge_dups

<details markdown="1">
<summary>Output files</summary>

- `purgedups/`
  - `<sample_id>.split.fasta`, `<sample_id>.split.self.paf.gz`: self-alignment inputs/outputs.
  - `<sample_id>.cutoffs`, `<sample_id>.hist.png`, `<sample_id>.calcuts.log`: `calcuts` thresholds and the coverage histogram plot for visual sanity-checking (not just the numbers).
  - `<sample_id>.dups.bed`, `<sample_id>.purge_dups.log`: detected duplication regions.
  - `<sample_id>.purged.fasta`: final purge_dups assembly. If `calcuts` could not determine usable cutoffs from this sample's coverage (not bimodal enough), this is the filtered assembly carried through unchanged, and the reason is recorded as a caveat in that sample's report.
  - `<sample_id>.purged_haplotigs.fasta`: contigs/regions removed as haplotigs.

</details>

### purge_haplotigs

<details markdown="1">
<summary>Output files</summary>

- `purgehaplotigs/`
  - `<sample_id>.bam.200.gencov`, `<sample_id>.bam.histogram.200.png`: `purge_haplotigs hist` output (`200` is the tool's `-d/-depth` cutoff, embedded in its output filenames).
  - `<sample_id>.cutoffs.json`, `<sample_id>.depth_hist.tsv`: automatically-estimated low/mid/high cutoffs and the depth histogram they were derived from. If no usable bimodal signal was found, `cutoffs.json` has `"skipped": true` and a `"reason"` instead of cutoffs, and the remaining files below are not produced for that sample -- the skip is instead recorded as a caveat in that sample's report.
  - `<sample_id>.coverage_stats.csv`: `purge_haplotigs cov` output.
  - `<sample_id>.curated.fasta`: final purge_haplotigs assembly.
  - `<sample_id>.curated.haplotigs.fasta`: contigs removed as haplotigs.

Only produced when `--run_purge_haplotigs` is `true` (default) and a usable bimodal coverage signal was found for that sample.

</details>

### Assembly stats

<details markdown="1">
<summary>Output files</summary>

- `assembly/<sample_id>.<stage>.stats.json`: contig count, total span, N50, longest contig, for `stage` in `raw`, `filtered`, `purge_dups`, `purge_haplotigs`.
- `assembly/<sample_id>.<stage>.assembly_stats_mqc.json`: the same numbers, reformatted as a MultiQC custom-content table (one row per sample x stage).

</details>

### BUSCO comparison

<details markdown="1">
<summary>Output files</summary>

- `filter/<sample_id>.raw.length_filtered.fasta`: the raw-stage assembly with contigs shorter than `--busco_min_contig_length` removed, used only as BUSCO's raw-stage input (not used anywhere else -- `assembly/` and the report's span/N50 always reflect the full, unfiltered raw assembly).
- `filter/<sample_id>.raw.contig_length_filter.json`, `<sample_id>.raw.contig_length_filter_mqc.json`: how many contigs/bp were excluded from the raw-stage BUSCO input and why.
- `split/<sample_id>.<stage>.chunk<K>of<N>.fasta`: each stage's (post-length-filter, for raw) assembly split into `N` span-balanced chunks per `--busco_chunk_max_span`, so no single BUSCO invocation has to hold the whole assembly's candidate alignments in memory at once.
- `split/<sample_id>.<stage>.split_summary.json`: how many chunks were produced and each chunk's contig count/span.
- `busco/<sample_id>.<stage>.<chunk_label>.<lineage>_busco/`: full BUSCO output directory per stage x chunk x lineage (`<chunk_label>` is e.g. `1of4`; always `1of1` when a stage wasn't split).
- `busco/<sample_id>.<stage>.<chunk_label>.<lineage>.short_summary.json`, `.full_table.tsv`: BUSCO's own per-chunk summary/gene-call table -- not meaningful on its own (each chunk only ever sees its own subset of contigs), kept for per-chunk debuggability.
- `merge/<sample_id>.<stage>.<lineage>.short_summary.json`: the per-ortholog merge of that (stage, lineage)'s chunk `full_table.tsv` files back into one Complete/Duplicated/Fragmented/Missing score for the whole stage -- this, not the raw per-chunk BUSCO output, is what feeds the comparison table and report below.
- `merge/short_summary.<sample_id>.<stage>.<lineage>.txt`: the same merged result in BUSCO's own summary-table text format, for MultiQC.
- `compare/<sample_id>.busco_comparison.tsv`, `<sample_id>.busco_comparison.json`: one comparative table per sample, across all stages and lineages, plus the per-lineage duplication-drop figures used in the report's verdict.
- `compare/<sample_id>.busco_comparison_mqc.json`: the same per-stage x lineage rows, reformatted as a MultiQC custom-content table.

</details>

Every requested `--busco_lineages` entry is run explicitly (never `--auto-lineage`) against the raw, filtered, purge_dups and (if enabled) purge_haplotigs assemblies -- except that BUSCO's raw-stage input is first length-filtered per `--busco_min_contig_length` (default 1000bp), and every stage's (post-filter) assembly is then split into span-balanced chunks per `--busco_chunk_max_span` (default 20Mbp) with results merged back per ortholog, both to avoid out-of-memory failures -- BUSCO's single-threaded post-processing scales with total genome content searched, not with fragmentation or duplication level, so even a clean, filtered assembly can exceed a memory ceiling if the chunk size isn't bounded. `assembly/` and the report's span/N50 figures are unaffected and always reflect the full, unfiltered raw assembly.

### blobpurge report

<details markdown="1">
<summary>Output files</summary>

- `analyses/<sample_id>_blobpurge.html`: the final, self-contained per-sample report -- span per stage, the BTK_FILTER span check, the GenomeScope2 comparison (or an explicit "not provided" notice), the BUSCO duplication trend per lineage, the purge_dups vs purge_haplotigs comparison (flagged if they disagree beyond `--purge_disagreement_threshold`), the collected caveats (ngscstat/Illumina path, purge_haplotigs on short-read coverage), and an explicit yes/partial/no verdict on whether uncollapsed heterozygous haplotigs explain the assembly's size/duplication surplus.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: aggregate report across all samples in the run (BUSCO results, samtools coverage stats, the BTK_FILTER span check, per-stage assembly stats, the cross-stage BUSCO comparison, and software versions).
  - `multiqc_data/`, `multiqc_plots/`: supporting data and static plot images.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`, `pipeline_dag.html`.
  - `nf-blobpurge_software_mqc_versions.yml`: software versions actually used, collated from every process.
  - `params_*.json`: parameters used by the pipeline run.

</details>
