# cnr-ibba/nf-blobpurge: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of cnr-ibba/nf-blobpurge, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- HaploMerger2 as a second, independent, opt-in haplotig-purging cross-check (`--run_haplomerger2`, default `false`), running in parallel to `purge_dups`/`purge_haplotigs` on the `BTK_FILTER` output. Unlike the other two methods, it purges by whole-genome self-alignment (LASTZ + UCSC kentUtils) rather than read coverage, so it needs no coverage/BAM input. New params: `--run_haplomerger2`, `--haplomerger2_identity`, `--haplomerger2_path`, `--haplomerger2_container` (mandatory when enabled -- no bioconda/biocontainers image exists for HaploMerger2, see `docs/usage.md`). The purge-method disagreement check in the per-sample report is now pairwise across however many purge methods actually ran (2 or 3), not hardcoded to purge_dups/purge_haplotigs.
- `docker/haplomerger2/Dockerfile`, a buildable recipe (verified against a real image built with `wave-cli`) providing HaploMerger2's external dependencies (`lastz` + 18 UCSC kentUtils tools + `perl`), since none of these ship together in any existing bioconda/biocontainers image.
- `HAPLOMERGER2_STAGE_A`/`HAPLOMERGER2_STAGE_B` now fail loudly with a pointed error if HaploMerger2 produces no usable sequence, instead of silently emitting an empty/corrupt assembly that would otherwise fail confusingly several steps downstream (e.g. in BUSCO). Verified to trigger correctly: HaploMerger2's `HM_pathFinder_preparation.pl`/`HM_pathFinder.pl` cannot build a workable alignment graph on assemblies with too few scaffolds/too little self-alignment signal, which is why `--run_haplomerger2` is not exercised by the default `-profile test` fixture yet (see `conf/test.config`).

### `Fixed`

### `Dependencies`

### `Deprecated`
