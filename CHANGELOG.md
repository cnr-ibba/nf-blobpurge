# cnr-ibba/nf-blobpurge: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of cnr-ibba/nf-blobpurge, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Support gzip-compressed per-field BlobDir JSON files (`*.json.gz`, as produced by newer BlobToolKit versions) in `bin/detect_taxrule.py` and `bin/check_span.py`.

### `Changed`

- Renamed `--run_organelle_isolation` to `--skip_organelle_isolation` and flipped its default so organelle isolation now runs **by default** (was opt-in): with it off, `BTK_FILTER` could misclassify and discard genuine mitochondrial/plastid contigs as bacterial contamination before organelle isolation ever ran, silently losing real organelle sequence from the delivered assembly. Set `--skip_organelle_isolation` to restore the old (isolation-off) behaviour.

### `Fixed`

- `purge_dups` and `purge_haplotigs` no longer abort the whole run when a sample's coverage isn't usefully bimodal (needed to auto-derive cutoffs): that sample's affected stage is now skipped individually, with a caveat explaining why in its report, instead of crashing the entire multi-sample pipeline.

### `Dependencies`

### `Deprecated`
