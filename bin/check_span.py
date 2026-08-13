#!/usr/bin/env python3
"""Verify that BTK_FILTER conserves assembly span.

Checks programmatically that:

    span(filtered assembly) + span(excluded contigs, from the BlobDir) == span(original assembly)

within a relative tolerance. The original and filtered spans are computed
directly from the FASTA files (ground truth of what was actually produced),
and the excluded-contig span is computed from the BlobDir's own length and
category fields (ground truth of what BTK_FILTER was asked to remove) -- not
from hand-entered numbers or by trusting a single tool's self-reported log.

Exits non-zero with an explicit message on mismatch, and writes a JSON report
either way so the discrepancy (or the pass) is traceable to numbers, not to
a bare pass/fail flag.
"""
import argparse
import gzip
import json
import sys
from pathlib import Path


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def resolve_blobdir_file(blobdir, filename):
    """Resolve <blobdir>/<filename>, preferring the gzip sibling
    <filename>.gz when both exist (newer BlobToolKit versions
    gzip-compress per-field BlobDir JSON files)."""
    base = Path(blobdir) / filename
    gz = Path(f"{base}.gz")
    if gz.exists():
        return gz
    if base.exists():
        return base
    raise FileNotFoundError(f"neither {base} nor {gz} exists")


def load_json(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def fasta_spans(path):
    """Return {record_id: length} for a (optionally gzipped) FASTA file."""
    spans = {}
    current_id = None
    current_len = 0
    with open_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    spans[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id is not None:
            spans[current_id] = current_len
    return spans


def field_values(data, path):
    """Extract a BlobDir field's value list from its parsed JSON dict."""
    for key in ("values", "data"):
        if key in data:
            return data[key]
    raise ValueError(f"{path}: no 'values' or 'data' key found in field JSON")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-fasta", required=True)
    parser.add_argument("--filtered-fasta", required=True)
    parser.add_argument("--blobdir", required=True, help="Path to the original (unfiltered) BlobDir")
    parser.add_argument("--taxon-field", required=True, help="Category field name, e.g. buscoregions_phylum")
    parser.add_argument("--exclude-taxa", required=True, help="Comma-separated taxa excluded by BTK_FILTER")
    parser.add_argument("--tolerance", type=float, default=0.001, help="Max relative discrepancy allowed")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument(
        "--exclude-ids",
        default=None,
        help=(
            "Optional newline-delimited contig ID file (e.g. bin/detect_organelles.py's "
            "--output-organelle-ids) to skip entirely when walking the BlobDir -- these "
            "contigs were isolated upstream of BTK_FILTER and are no longer part of "
            "--original-fasta's scope, so they must not be double-counted into excluded_span."
        ),
    )
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-mqc-json", required=True, help="MultiQC custom-content JSON path")
    args = parser.parse_args()

    exclude_taxa = {t.strip() for t in args.exclude_taxa.split(",") if t.strip()}

    skip_ids = set()
    if args.exclude_ids:
        with open(args.exclude_ids) as handle:
            skip_ids = {line.strip() for line in handle if line.strip()}

    original_spans = fasta_spans(args.original_fasta)
    filtered_spans = fasta_spans(args.filtered_fasta)
    original_span = sum(original_spans.values())
    filtered_span = sum(filtered_spans.values())

    identifiers = field_values(load_json(resolve_blobdir_file(args.blobdir, "identifiers.json")), "identifiers.json")
    lengths = field_values(load_json(resolve_blobdir_file(args.blobdir, "length.json")), "length.json")

    try:
        taxon_path = resolve_blobdir_file(args.blobdir, f"{args.taxon_field}.json")
    except FileNotFoundError:
        sys.exit(
            f"ERROR: BlobDir field '{args.taxon_field}.json' (or '.json.gz') not found in {args.blobdir}. "
            "--taxon_field must name a categorical field that actually exists in this "
            "BlobDir (check <blobdir>/meta.json for the list of available 'fields')."
        )
    field_meta = load_json(taxon_path)
    categories = field_values(field_meta, taxon_path)
    keys = field_meta.get("keys")

    if len(identifiers) != len(lengths) or len(identifiers) != len(categories):
        sys.exit(
            "ERROR: BlobDir field arrays have inconsistent lengths "
            f"(identifiers={len(identifiers)}, length={len(lengths)}, "
            f"{args.taxon_field}={len(categories)}). This BlobDir looks corrupt "
            "or was not built from the same assembly given in the samplesheet."
        )

    excluded_span = 0
    excluded_ids = []
    retained_ids = []
    for ctg_id, ctg_len, cat in zip(identifiers, lengths, categories):
        if ctg_id in skip_ids:
            continue
        label = keys[cat] if keys is not None and isinstance(cat, int) else cat
        if label in exclude_taxa:
            excluded_span += ctg_len
            excluded_ids.append(ctg_id)
        else:
            retained_ids.append(ctg_id)

    reconstructed = filtered_span + excluded_span
    diff = abs(reconstructed - original_span)
    relative_diff = diff / original_span if original_span else 0.0
    passed = relative_diff <= args.tolerance

    report = {
        "sample_id": args.sample_id,
        "original_span": original_span,
        "original_n_contigs": len(original_spans),
        "filtered_span": filtered_span,
        "filtered_n_contigs": len(filtered_spans),
        "excluded_span_from_blobdir": excluded_span,
        "excluded_n_contigs": len(excluded_ids),
        "reconstructed_span": reconstructed,
        "absolute_diff": diff,
        "relative_diff": relative_diff,
        "tolerance": args.tolerance,
        "taxon_field": args.taxon_field,
        "exclude_taxa": sorted(exclude_taxa),
        "excluded_ids_skipped": sorted(skip_ids),
        "excluded_ids_skipped_n_contigs": len(skip_ids),
        "pass": passed,
    }

    with open(args.output_json, "w") as handle:
        json.dump(report, handle, indent=2)

    mqc = {
        "id": "btk_filter_span_check",
        "section_name": "Contamination filtering (span check)",
        "description": (
            "Programmatic check that span(filtered) + span(excluded, from the BlobDir) "
            "equals span(original), within --span_tolerance."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "btk_filter_span_check_table",
            "title": "BTK_FILTER span check",
            "namespace": "BTK Filter",
        },
        "data": {
            args.sample_id: {
                "Sample": args.sample_id,
                "Pass": "Yes" if passed else "No",
                "Original span (bp)": original_span,
                "Filtered span (bp)": filtered_span,
                "Excluded span (bp)": excluded_span,
                "Relative diff": relative_diff,
            }
        },
    }
    with open(args.output_mqc_json, "w") as handle:
        json.dump(mqc, handle, indent=2)

    with open("retained_ids.txt", "w") as handle:
        handle.write("\n".join(retained_ids) + ("\n" if retained_ids else ""))
    with open("excluded_ids.txt", "w") as handle:
        handle.write("\n".join(excluded_ids) + ("\n" if excluded_ids else ""))

    if not passed:
        sys.exit(
            "ERROR: BTK_FILTER span check failed for sample "
            f"'{args.sample_id}': filtered assembly span ({filtered_span}) + "
            f"excluded contigs span from BlobDir ({excluded_span}) = "
            f"{reconstructed}, vs original assembly span {original_span} "
            f"(relative diff {relative_diff:.4%} > tolerance {args.tolerance:.4%}). "
            "This usually means --taxon_field/--taxrule do not match how the "
            "BlobDir was built, or the assembly given in the samplesheet is not "
            "the one the BlobDir was generated from. See docs/usage.md."
        )

    print(f"BTK_FILTER span check passed for '{args.sample_id}' (relative diff {relative_diff:.6%}).")


if __name__ == "__main__":
    main()
