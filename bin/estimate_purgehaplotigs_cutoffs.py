#!/usr/bin/env python3
"""Auto-estimate low/mid/high read-depth cutoffs for `purge_haplotigs cov`.

purge_haplotigs normally expects a human to eyeball the coverage histogram
produced by `purge_haplotigs hist` and pick three cutoffs by hand. That is
not compatible with a generic, non-interactive pipeline, so this script
derives them from the genome-wide depth distribution instead (from
`samtools depth -a`), assuming the classic bimodal signature of an assembly
with uncollapsed heterozygous haplotigs: a "haploid" peak at roughly half
the depth of the main "diploid/collapsed" peak.

This is a heuristic over the sample's own data (not a per-species constant):
    - primary peak  = the modal depth (diploid/collapsed coverage)
    - haploid peak  = the modal depth in the [0.25, 0.75] * primary window
    - low           = the depth-histogram trough before the haploid peak
    - mid           = the depth-histogram trough between the two peaks
    - high          = primary peak * --high-multiplier

If no distinct haploid peak is found (e.g. an already well-purged or haploid
assembly), the script does not guess: it writes a JSON result with
"skipped": true and an explanation instead of low/mid/high cutoffs, and exits
0 -- the caller (the Nextflow subworkflow) skips purge_haplotigs cov/purge for
this sample rather than treating "no bimodal signal" as a pipeline failure.
"""
import argparse
import json
import sys
from collections import Counter


def read_depth_histogram(depth_tsv, max_depth_cap):
    counts = Counter()
    with open(depth_tsv) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            depth = int(fields[2])
            if depth <= 0:
                continue
            counts[min(depth, max_depth_cap)] += 1
    return counts


def smooth(counts, max_depth, window=3):
    smoothed = {}
    for d in range(1, max_depth + 1):
        lo, hi = max(1, d - window), min(max_depth, d + window)
        vals = [counts.get(x, 0) for x in range(lo, hi + 1)]
        smoothed[d] = sum(vals) / len(vals)
    return smoothed


def find_peak(smoothed, lo, hi):
    window = {d: v for d, v in smoothed.items() if lo <= d <= hi}
    if not window:
        return None
    return max(window, key=window.get)


def find_trough(smoothed, lo, hi):
    window = {d: v for d, v in smoothed.items() if lo <= d <= hi}
    if not window:
        return None
    return min(window, key=window.get)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("depth_tsv", help="Output of `samtools depth -a aln.bam`")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--max-depth-cap", type=int, default=500, help="Depths above this are folded into the cap bin (repeat-region tail)")
    parser.add_argument("--high-multiplier", type=float, default=2.0)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--histogram-tsv", required=True)
    args = parser.parse_args()

    def write_skipped(reason, **extra):
        report = {"sample_id": args.sample_id, "skipped": True, "reason": reason, **extra}
        with open(args.output_json, "w") as handle:
            json.dump(report, handle, indent=2)
        print(f"purge_haplotigs cutoffs skipped for '{args.sample_id}': {reason}", file=sys.stderr)
        print(json.dumps(report, indent=2))

    counts = read_depth_histogram(args.depth_tsv, args.max_depth_cap)
    if not counts:
        open(args.histogram_tsv, "w").close()
        write_skipped(f"no non-zero depth positions found in {args.depth_tsv}; cannot estimate purge_haplotigs cutoffs.")
        return

    max_depth = max(counts)
    smoothed = smooth(counts, max_depth)

    with open(args.histogram_tsv, "w") as handle:
        handle.write("depth\tn_positions\tsmoothed\n")
        for d in range(1, max_depth + 1):
            handle.write(f"{d}\t{counts.get(d, 0)}\t{smoothed[d]:.2f}\n")

    primary_peak = find_peak(smoothed, 1, max_depth)
    haploid_lo = max(1, int(primary_peak * 0.25))
    haploid_hi = max(haploid_lo + 1, int(primary_peak * 0.75))
    haploid_peak = find_peak(smoothed, haploid_lo, haploid_hi)

    if haploid_peak is None or haploid_peak >= primary_peak:
        write_skipped(
            f"could not identify a distinct haploid coverage peak below the primary peak "
            f"({primary_peak}x). The coverage distribution does not show the classic bimodal "
            "signature of uncollapsed heterozygous haplotigs; purge_haplotigs cutoffs cannot "
            f"be estimated automatically. Inspect {args.histogram_tsv} and, if appropriate, "
            "supply cutoffs manually.",
            primary_peak=primary_peak,
        )
        return

    low = find_trough(smoothed, max(1, int(haploid_peak * 0.3)), haploid_peak)
    mid = find_trough(smoothed, haploid_peak, primary_peak)
    high = min(max_depth, int(primary_peak * args.high_multiplier))

    if low is None or mid is None or not (low < haploid_peak < mid < primary_peak < high):
        write_skipped(
            f"depth troughs around the candidate haploid peak ({haploid_peak}x) and primary "
            f"peak ({primary_peak}x) were not well separated; refusing to guess purge_haplotigs "
            f"cutoffs. Inspect {args.histogram_tsv} and supply cutoffs manually.",
            primary_peak=primary_peak,
            haploid_peak=haploid_peak,
        )
        return

    report = {
        "sample_id": args.sample_id,
        "skipped": False,
        "primary_peak": primary_peak,
        "haploid_peak": haploid_peak,
        "low": low,
        "mid": mid,
        "high": high,
    }
    with open(args.output_json, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
