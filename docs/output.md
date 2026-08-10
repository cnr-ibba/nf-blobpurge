# cnr-ibba/nf-blobpurge: Output

## Introduction

This document describes the output produced by the pipeline for each sample in the samplesheet. All paths below are relative to the top-level `--outdir` results directory.

## Pipeline overview

The pipeline processes each sample through the following stages:

- [BTK_FILTER](#btk_filter) - contamination removal from an existing BlobDir, with a programmatic span-conservation check
- [Read coverage](#read-coverage) - CRAM subsetting or bwa-mem2 mapping, feeding purge_dups' `ngscstat`
- [Organelle isolation](#organelle-isolation) - optional coverage-based isolation of mitochondrial/plastid contigs before purging
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
- `samtools/<sample_id>.namesorted.bam`: the same BAM, name-sorted, feeding `ngscstat` below.
- `samtools/<sample_id>.nuclear.namesorted.bam`, `<sample_id>.organelle.namesorted.bam`: only present with `--run_organelle_isolation`; the nuclear-only and organelle-only BAM subsets, name-sorted, feeding the organelle-free `ngscstat` re-derivation and organelle read export (see [Organelle isolation](#organelle-isolation)).
- `purgedups/<sample_id>.ngscstat.stat`, `<sample_id>.ngscstat.base.cov`: purge_dups Illumina coverage stats (`ngscstat`).

</details>

### Organelle isolation

<details markdown="1">
<summary>Output files</summary>

- `organelle/`
  - `<sample_id>.organelle_report.json`: classification report -- baseline/threshold coverage used, and, per contig, mean depth, coverage ratio, length, and (if `--organelle_reference_fasta` was given) reference-match fraction.
  - `<sample_id>.organelle_report.mqc.json`: the same summary, reformatted as a MultiQC custom-content table.
  - `<sample_id>.organelle.fasta` / `<sample_id>.nuclear.fasta`: the filtered assembly split by classification.
  - `<sample_id>.organelle_ids.txt` / `<sample_id>.nuclear_ids.txt`: the corresponding contig ID lists.
  - `<sample_id>.organelle_ref.paf`: minimap2 alignment of the filtered assembly against `--organelle_reference_fasta`, if given (empty otherwise).
  - `<sample_id>.organelle_R1.fastq.gz` / `<sample_id>.organelle_R2.fastq.gz`: paired-end reads mapping to the isolated organelle contigs, for use with external organelle-assembly tools (this pipeline does not assemble organelles itself).
  - `<sample_id>.purge_dups.merged.fasta` / `<sample_id>.purge_haplotigs.merged.fasta`: the isolated organelle contigs re-merged back into each purge stage's assembly.
- `nuclear/`: the coordinate-sorted/indexed BAM subset to nuclear-only contigs, used to re-derive an organelle-free `ngscstat` coverage pair (see [purge_dups](#purge_dups)) and fed to `purge_haplotigs`.

</details>

Only produced when `--run_organelle_isolation` is `true` (default `false`). Isolated contigs are not contamination and are not discarded: `purge_dups`/`purge_haplotigs` run on the nuclear-only assembly so their coverage-cutoff estimation isn't distorted by the organelle's much higher copy number, and the isolated contigs are re-merged into the final assembly at every stage afterwards.

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

- `busco/<sample_id>.<stage>.<lineage>_busco/`: full BUSCO output directory per stage x lineage.
- `busco/<sample_id>.<stage>.<lineage>.short_summary.json`: BUSCO's own summary JSON, the source of every number used downstream.
- `filter/<sample_id>.raw.length_filtered.fasta`: the raw-stage assembly with contigs shorter than `--busco_min_contig_length` removed, used only as BUSCO's raw-stage input (not used anywhere else -- `assembly/` and the report's span/N50 always reflect the full, unfiltered raw assembly).
- `filter/<sample_id>.raw.contig_length_filter.json`, `<sample_id>.raw.contig_length_filter_mqc.json`: how many contigs/bp were excluded from the raw-stage BUSCO input and why.
- `compare/<sample_id>.busco_comparison.tsv`, `<sample_id>.busco_comparison.json`: one comparative table per sample, across all stages and lineages, plus the per-lineage duplication-drop figures used in the report's verdict.
- `compare/<sample_id>.busco_comparison_mqc.json`: the same per-stage x lineage rows, reformatted as a MultiQC custom-content table.

</details>

Every requested `--busco_lineages` entry is run explicitly (never `--auto-lineage`) against the raw, filtered, purge_dups and (if enabled) purge_haplotigs assemblies -- except that BUSCO's raw-stage input is first length-filtered per `--busco_min_contig_length` (default 1000bp) to avoid an out-of-memory failure on highly fragmented assemblies; `assembly/` and the report's span/N50 figures are unaffected and always reflect the full, unfiltered raw assembly.

### blobpurge report

<details markdown="1">
<summary>Output files</summary>

- `analyses/<sample_id>_blobpurge.html`: the final, self-contained per-sample report -- span per stage, the BTK_FILTER span check, organelle contig isolation (or an explicit "not enabled" notice), the GenomeScope2 comparison (or an explicit "not provided" notice), the BUSCO duplication trend per lineage, the purge_dups vs purge_haplotigs comparison (flagged if they disagree beyond `--purge_disagreement_threshold`), the collected caveats (ngscstat/Illumina path, purge_haplotigs on short-read coverage), and an explicit yes/partial/no verdict on whether uncollapsed heterozygous haplotigs explain the assembly's size/duplication surplus.

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
