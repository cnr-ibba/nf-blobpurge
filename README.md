# cnr-ibba/nf-blobpurge

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/cnr-ibba/nf-blobpurge)
[![GitHub Actions CI Status](https://github.com/cnr-ibba/nf-blobpurge/actions/workflows/nf-test.yml/badge.svg)](https://github.com/cnr-ibba/nf-blobpurge/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/cnr-ibba/nf-blobpurge/actions/workflows/linting.yml/badge.svg)](https://github.com/cnr-ibba/nf-blobpurge/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/cnr-ibba/nf-blobpurge)

## Introduction

**cnr-ibba/nf-blobpurge** is a bioinformatics pipeline for de novo short-read genome assemblies (e.g. algae) that have already been classified with [sanger-tol/blobtoolkit](https://github.com/sanger-tol/blobtoolkit). Starting from an existing assembly + BlobDir pair (it never regenerates the BlobDir), it removes contaminant contigs, purges uncollapsed heterozygous haplotigs with `purge_dups` (cross-checked independently with `purge_haplotigs`), runs comparative BUSCO across every assembly stage, and produces a single per-sample HTML report with an explicit verdict on whether the size/duplication surplus is explained by heterozygous haplotigs.

Nothing about a specific organism, taxon, or BUSCO lineage is hardcoded: the taxa to exclude and the BUSCO lineage(s) to run are mandatory parameters with no default.

1. Remove contaminant contigs from an existing BlobDir, with a programmatic assembly-span conservation check (`BTK_FILTER`)
2. Compute read coverage for purging: subset the existing reads CRAM produced upstream by blobtoolkit (`READ_COVERAGE`)
3. Purge uncollapsed heterozygous haplotigs (`purge_dups`), cross-checked independently and in parallel with `purge_haplotigs`
4. Run comparative BUSCO (never `--auto-lineage`) across the raw, filtered and purged assemblies
5. Generate a per-sample HTML report tying span, BUSCO duplication, the GenomeScope2 comparison (if provided) and the purge_dups/purge_haplotigs cross-check together into an explicit verdict
6. Present QC and software versions for the whole run ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows (see [docs/usage.md](docs/usage.md) for the full column reference):

`samplesheet.csv`:

```csv
sample_id,assembly,blobdir,reads_cram
sample1,/data/sample1.assembly.fasta,/data/sample1_blobdir,/data/sample1.reads.cram
```

Each row represents one already-assembled, already-blobtoolkit-classified sample: its assembly FASTA, its BlobDir, and the reads CRAM blobtoolkit already produced (aligned to the raw assembly) for coverage.

Now, you can run the pipeline using:

```bash
nextflow run cnr-ibba/nf-blobpurge \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --exclude_taxa "Pseudomonadota,Bacteroidota,Actinomycetota" \
   --busco_lineages "chlorophyta_odb12,viridiplantae_odb12"
```

`--exclude_taxa` and `--busco_lineages` are mandatory and have no default -- the pipeline fails immediately with an explicit error if either is missing.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

cnr-ibba/nf-blobpurge was originally written by Paolo Cozzi.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use cnr-ibba/nf-blobpurge for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
