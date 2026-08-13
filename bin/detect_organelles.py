#!/usr/bin/env python3
"""Detect organelle-like (mitochondrial/plastid) contigs by extreme coverage.

Organelle genomes have a much higher copy number per cell than the nuclear
genome, so their sequencing depth is typically many times the nuclear
baseline. Left in the assembly fed to PURGE_DUPS/PURGE_HAPLOTIGS, that extreme
coverage can distort the genome-wide coverage histograms those tools use to
auto-estimate cutoffs (both `calcuts` and `estimate_purgehaplotigs_cutoffs.py`
pool coverage across all contigs, with no per-contig awareness).

This script classifies contigs from a per-contig coverage table (the output
of `samtools coverage`) rather than trying to fix the pooled histograms after
the fact:

    baseline  = length-weighted median of mean depth across all contigs
    threshold = --coverage-multiplier * baseline
    isolated  = (mean_depth >= threshold) AND (length <= --max-length, if set)

A reference FASTA match (via a minimap2 PAF against known mitochondrial/
plastid sequences) is optional corroborating evidence, recorded on every
contig regardless of the coverage-based decision. By default it is purely
informational; pass --require-reference-match to make it a hard requirement
in addition to the coverage/length gate.

Exits 0 in every case, including zero contigs isolated -- that is a valid
biological outcome (e.g. no reads recruited to an organelle, or an organelle
already collapsed with the nuclear assembly), not a failure.
"""
import argparse
import gzip
import json
import sys
from collections import defaultdict


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_coverage_table(path):
    """Parse `samtools coverage` output: #rname startpos endpos numreads
    covbases coverage meandepth meanbaseq meanmapq (tab-separated)."""
    contigs = {}
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            rname, endpos, meandepth = fields[0], int(fields[2]), float(fields[6])
            contigs[rname] = {"length": endpos, "mean_depth": meandepth}
    return contigs


def parse_paf_aligned_fractions(path):
    """Return {query_id: aligned_fraction} from a minimap2 PAF, merging
    overlapping query-side alignment intervals so a contig hit by several
    overlapping alignments isn't counted more than once."""
    if path is None:
        return {}
    intervals = defaultdict(list)
    qlengths = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            qlengths[qname] = qlen
            intervals[qname].append((qstart, qend))

    fractions = {}
    for qname, ivs in intervals.items():
        ivs.sort()
        merged = []
        for start, end in ivs:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        aligned = sum(e - s for s, e in merged)
        qlen = qlengths[qname]
        fractions[qname] = aligned / qlen if qlen else 0.0
    return fractions


def weighted_median(values_weights):
    """Length-weighted median: robust to the very high-coverage outliers
    this script is trying to detect (a plain mean would be dragged up by
    them; the median of a handful of contigs would be dominated by whichever
    contig happens to be picked, regardless of how much sequence it covers)."""
    items = sorted(v for v in values_weights if v[1] > 0)
    total = sum(w for _, w in items)
    if total == 0:
        return 0.0
    half = total / 2.0
    cum = 0.0
    for value, weight in items:
        cum += weight
        if cum >= half:
            return value
    return items[-1][0]


def split_fasta(assembly_fasta, organelle_ids, organelle_out, nuclear_out):
    with open_maybe_gzip(assembly_fasta) as handle, \
            open(organelle_out, "w") as org_handle, \
            open(nuclear_out, "w") as nuc_handle:
        out = nuc_handle
        for line in handle:
            if line.startswith(">"):
                contig_id = line[1:].split()[0].strip()
                out = org_handle if contig_id in organelle_ids else nuc_handle
            out.write(line if line.endswith("\n") else line + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--assembly-fasta", required=True)
    parser.add_argument("--coverage-table", required=True, help="Output of `samtools coverage`")
    parser.add_argument("--reference-paf", default=None, help="Optional minimap2 PAF of assembly contigs vs an organelle reference")
    parser.add_argument("--coverage-multiplier", type=float, required=True)
    parser.add_argument("--max-length", type=int, default=None, help="Optional AND-gate: only contigs at or below this length are eligible")
    parser.add_argument("--min-reference-coverage", type=float, default=0.8)
    parser.add_argument("--require-reference-match", action="store_true", help="Require a reference match in addition to the coverage/length gate")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-mqc-json", required=True)
    parser.add_argument("--output-caveat", required=True)
    parser.add_argument("--output-organelle-ids", required=True)
    parser.add_argument("--output-nuclear-ids", required=True)
    parser.add_argument("--output-organelle-fasta", required=True)
    parser.add_argument("--output-nuclear-fasta", required=True)
    args = parser.parse_args()

    contigs = parse_coverage_table(args.coverage_table)
    if not contigs:
        sys.exit(f"ERROR: no contigs found in coverage table '{args.coverage_table}'.")

    reference_used = args.reference_paf is not None
    reference_fractions = parse_paf_aligned_fractions(args.reference_paf)

    baseline = weighted_median([(c["mean_depth"], c["length"]) for c in contigs.values()])
    threshold = args.coverage_multiplier * baseline
    # A zero baseline (e.g. reads recruit to none of the contigs' majority
    # span -- plausible here since this runs on the RAW, unfiltered assembly,
    # before BTK_FILTER has removed anything the reads don't map to) makes
    # threshold 0 too, so `mean_depth >= threshold` would be true for every
    # contig and the whole assembly would be misclassified as organelle.
    # Treat "no usable coverage signal" as "isolate nothing" instead.
    no_coverage_signal = baseline <= 0

    rows = []
    organelle_ids = []
    nuclear_ids = []
    total_isolated_span = 0
    n_coverage_candidates = 0
    for contig_id in sorted(contigs):
        info = contigs[contig_id]
        length, mean_depth = info["length"], info["mean_depth"]
        coverage_flag = False if no_coverage_signal else mean_depth >= threshold
        length_flag = True if args.max_length is None else length <= args.max_length
        isolated = coverage_flag and length_flag
        if coverage_flag:
            n_coverage_candidates += 1

        reference_aligned_fraction = None
        reference_matched = None
        if reference_used:
            reference_aligned_fraction = reference_fractions.get(contig_id, 0.0)
            reference_matched = reference_aligned_fraction >= args.min_reference_coverage
            if args.require_reference_match:
                isolated = isolated and reference_matched

        rows.append({
            "id": contig_id,
            "length": length,
            "mean_depth": mean_depth,
            "coverage_ratio": (mean_depth / baseline) if baseline else None,
            "coverage_flag": coverage_flag,
            "length_flag": length_flag,
            "reference_aligned_fraction": reference_aligned_fraction,
            "reference_matched": reference_matched,
            "isolated": isolated,
        })

        if isolated:
            organelle_ids.append(contig_id)
            total_isolated_span += length
        else:
            nuclear_ids.append(contig_id)

    report = {
        "sample_id": args.sample_id,
        "enabled": True,
        "baseline_coverage": baseline,
        "coverage_multiplier": args.coverage_multiplier,
        "threshold": threshold,
        "no_coverage_signal": no_coverage_signal,
        "max_length": args.max_length,
        "reference_fasta_used": reference_used,
        "min_reference_coverage": args.min_reference_coverage if reference_used else None,
        "require_reference_match": args.require_reference_match,
        "n_contigs_total": len(contigs),
        "n_coverage_candidates": n_coverage_candidates,
        "n_isolated": len(organelle_ids),
        "total_isolated_span": total_isolated_span,
        "contigs": rows,
    }
    with open(args.output_json, "w") as handle:
        json.dump(report, handle, indent=2)

    mqc = {
        "id": "organelle_isolation",
        "section_name": "Organelle contig isolation",
        "description": (
            "Contigs classified as organelle-like by coverage "
            f"(>= {args.coverage_multiplier}x the length-weighted median contig depth)"
            + (", optionally corroborated by a reference match" if reference_used else "")
            + ", and isolated from purge_dups/purge_haplotigs coverage-cutoff estimation."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "organelle_isolation_table",
            "title": "Organelle contig isolation",
            "namespace": "Organelle isolation",
        },
        "data": {
            args.sample_id: {
                "Sample": args.sample_id,
                "Baseline coverage": round(baseline, 2),
                "Threshold": round(threshold, 2),
                "Contigs isolated": len(organelle_ids),
                "Isolated span (bp)": total_isolated_span,
            }
        },
    }
    with open(args.output_mqc_json, "w") as handle:
        json.dump(mqc, handle, indent=2)

    if no_coverage_signal:
        caveat = (
            f"Organelle isolation enabled for sample '{args.sample_id}' but the length-weighted "
            "median contig depth was 0x (no usable coverage signal, e.g. reads did not map to a "
            "majority of the assembly's span): no contigs were isolated rather than risk "
            "misclassifying the whole assembly as organelle-like."
        )
    elif organelle_ids:
        caveat = (
            f"{len(organelle_ids)} organelle-like contig(s) isolated for sample '{args.sample_id}' "
            f"(total {total_isolated_span} bp, mean depth >= {threshold:.1f}x, "
            f"{args.coverage_multiplier}x the {baseline:.1f}x nuclear baseline): "
            f"{', '.join(organelle_ids)}. Excluded from purge_dups/purge_haplotigs coverage-cutoff "
            "estimation and re-merged into the final assembly at each stage; a matching paired-end "
            "read subset was exported for use with external organelle-assembly tools "
            "(e.g. GetOrganelle, MitoHiFi, oatk)."
        )
    elif n_coverage_candidates:
        caveat = (
            f"Organelle isolation enabled for sample '{args.sample_id}': {n_coverage_candidates} "
            f"contig(s) reached the coverage threshold ({args.coverage_multiplier}x the "
            f"{baseline:.1f}x nuclear baseline) but were excluded by --organelle_max_length and/or "
            "--organelle_require_reference_match: no contigs were isolated."
        )
    else:
        caveat = (
            f"Organelle isolation enabled for sample '{args.sample_id}' but no contig reached "
            f"the coverage threshold ({args.coverage_multiplier}x the {baseline:.1f}x nuclear "
            "baseline): no contigs were isolated."
        )
    with open(args.output_caveat, "w") as handle:
        handle.write(caveat + "\n")

    with open(args.output_organelle_ids, "w") as handle:
        handle.write("\n".join(organelle_ids) + ("\n" if organelle_ids else ""))
    with open(args.output_nuclear_ids, "w") as handle:
        handle.write("\n".join(nuclear_ids) + ("\n" if nuclear_ids else ""))

    split_fasta(args.assembly_fasta, set(organelle_ids), args.output_organelle_fasta, args.output_nuclear_fasta)

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
