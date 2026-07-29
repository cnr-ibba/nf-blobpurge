# cnr-ibba/nf-blobpurge: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of cnr-ibba/nf-blobpurge, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- HaploMerger2 as a second, independent, opt-in haplotig-purging cross-check (`--run_haplomerger2`, default `false`), running in parallel to `purge_dups`/`purge_haplotigs` on the `BTK_FILTER` output. Unlike the other two methods, it purges by whole-genome self-alignment (LASTZ + UCSC kentUtils) rather than read coverage, so it needs no coverage/BAM input. New params: `--run_haplomerger2`, `--haplomerger2_identity`, `--haplomerger2_path`, `--haplomerger2_container` (mandatory when enabled -- no bioconda/biocontainers image exists for HaploMerger2, see `docs/usage.md`). The purge-method disagreement check in the per-sample report is now pairwise across however many purge methods actually ran (2 or 3), not hardcoded to purge_dups/purge_haplotigs.
- `docker/haplomerger2/Dockerfile`, a buildable recipe (verified against a real image built with `wave-cli`) providing HaploMerger2's external dependencies (`lastz` + 18 UCSC kentUtils tools + `perl`), since none of these ship together in any existing bioconda/biocontainers image.
- `HAPLOMERGER2_STAGE_A`/`HAPLOMERGER2_STAGE_B` now fail loudly with a pointed error if HaploMerger2 produces no usable sequence, instead of silently emitting an empty/corrupt assembly that would otherwise fail confusingly several steps downstream (e.g. in BUSCO).
- Grew the synthetic `-profile test` fixture from 6 to 9 contigs (a second near-duplicate pair, plus a larger unique host contig to preserve the single-copy-vs-haplotig depth ratio `purge_haplotigs`' cutoff estimation relies on) to exercise HaploMerger2's self-alignment graph more thoroughly, and switched `assets/test/assembly.fasta` to gzipped `assembly.fasta.gz` (the samplesheet schema already accepted `.fasta.gz`; reads were already gzipped).

### `Fixed`

- `BUSCO` (`modules/local/busco/main.nf`) now decompresses its input FASTA before running (`zcat -f`), since BUSCO/Biopython does not auto-detect gzip -- needed once the `raw` stage started receiving `assembly.fasta.gz` directly from the samplesheet.
- `HAPLOMERGER2_STAGE_A`/`HAPLOMERGER2_STAGE_B` no longer inherit Nextflow's own strict shell options (`errexit`/`nounset`/`pipefail`, auto-exported as `SHELLOPTS`) into HaploMerger2's `#!/bin/bash` `hm.batch*` scripts (via `env -u SHELLOPTS`) -- HaploMerger2 was never written to tolerate that strictness, and inheriting it caused silent mid-script aborts that only surfaced several steps later as a confusing "no usable sequence" failure.
- `HAPLOMERGER2_STAGE_A`/`HAPLOMERGER2_STAGE_B` no longer put the whole `chainNet_jksrc20100603_centOS5/` directory (needed only for the unpackaged `faToNib`) on `PATH` -- it also bundles ancient 2010 copies of `axtChain`/`chainNet`/etc. that were silently shadowing the modern, bioconda-installed versions from the container. Only `faToNib` itself is exposed now, via a dedicated directory.
- Known, unresolved: even with both fixes above, HaploMerger2's reciprocal-best chain/net computation (`hm.batchA2`) has been observed to intermittently produce empty output with no error, which then looks like a too-sparse alignment graph downstream. Root cause not isolated; `--run_haplomerger2` therefore remains `false` in `conf/test.config` and unexercised by CI -- see `docs/usage.md`.
- `tests/.nftignore` now excludes `busco/**` and `purgehaplotigs/*.bam.200.gencov` from the nf-test content snapshot: BUSCO's own `short_summary.json`/logs/per-marker hmmsearch output embed the absolute Nextflow task work-dir path (unique per run by construction, not a real BUSCO result difference), and `purge_haplotigs hist` writes its per-contig coverage blocks in non-reproducible order across runs (values identical, order only -- not a threading race, confirmed by pinning `-t 1`). Both were previously masked by the smaller 6-contig fixture and only became visible once it grew to 9 contigs.

### `Dependencies`

### `Deprecated`
