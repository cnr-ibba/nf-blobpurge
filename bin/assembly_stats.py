#!/usr/bin/env python3
"""Compute basic span/contig-count stats for a FASTA file, as JSON.

Used to record a traceable, machine-generated span number for each pipeline
stage (raw / filtered / purge_dups / purge_haplotigs), so the final report
never has to re-derive or hand-enter these figures.
"""
import argparse
import gzip
import json


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage", required=True, help="e.g. raw, filtered, purge_dups, purge_haplotigs")
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    lengths = []
    current_len = 0
    started = False
    with open_maybe_gzip(args.fasta) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if started:
                    lengths.append(current_len)
                current_len = 0
                started = True
            else:
                current_len += len(line.strip())
        if started:
            lengths.append(current_len)

    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = None
    running = 0
    for length in lengths:
        running += length
        if running >= total / 2:
            n50 = length
            break

    report = {
        "sample_id": args.sample_id,
        "stage": args.stage,
        "n_contigs": len(lengths),
        "total_span": total,
        "n50": n50,
        "longest_contig": lengths[0] if lengths else 0,
    }
    with open(args.output_json, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
