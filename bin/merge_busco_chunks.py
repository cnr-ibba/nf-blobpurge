#!/usr/bin/env python3
"""Merge per-chunk BUSCO full_table.tsv files into one ortholog-level result.

When an assembly stage is split into span-balanced chunks (see
bin/split_assembly_for_busco.py) and BUSCO is run independently on each
chunk, a true duplicated BUSCO ortholog may have its two copies land in two
different chunks -- each per-chunk BUSCO run only ever sees its own subset
of contigs, so it can only call a busco "Duplicated" on its own if both
copies happen to fall in the same chunk. This script pools every chunk's
full_table.tsv rows per busco_id and re-applies BUSCO's own classification
logic (count of Complete/Duplicated matches: 0 -> Missing/Fragmented,
1 -> Complete, >=2 -> Duplicated) across the whole assembly. This is exact,
not an approximation: chunks are always whole, non-overlapping sets of
contigs, so no individual gene call is ever split across chunks.

Writes a short_summary.json shaped exactly like BUSCO's own
short_summary...json (the subset of fields bin/compare_busco.py reads), and
a minimal MultiQC-parseable short_summary...txt, so this script's output is
a drop-in replacement for BUSCO's own per-run output in every downstream
consumer.
"""
import argparse
import re
from collections import defaultdict

HEADER_RE = re.compile(r"number of BUSCOs:\s*(\d+)")
FULL_MATCH_STATUSES = {"Complete", "Duplicated"}


def parse_full_table(path):
    """Yields (busco_id, status) for every data row, and the dataset's total
    BUSCO count read from the file's own header (constant across chunks)."""
    n_markers = None
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if n_markers is None:
                    match = HEADER_RE.search(line)
                    if match:
                        n_markers = int(match.group(1))
                continue
            fields = line.split("\t")
            rows.append((fields[0], fields[1]))
    return rows, n_markers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("full_tables", nargs="+", help="Per-chunk full_table.tsv files for one (sample, stage, lineage)")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-txt", required=True)
    args = parser.parse_args()

    statuses_by_busco = defaultdict(list)
    n_markers = None
    for path in args.full_tables:
        rows, header_n_markers = parse_full_table(path)
        if header_n_markers is not None:
            n_markers = header_n_markers
        for busco_id, status in rows:
            statuses_by_busco[busco_id].append(status)

    if n_markers is None:
        n_markers = len(statuses_by_busco)

    single_copy = 0
    multi_copy = 0
    fragmented = 0
    for statuses in statuses_by_busco.values():
        n_full = sum(1 for s in statuses if s in FULL_MATCH_STATUSES)
        if n_full == 1:
            single_copy += 1
        elif n_full >= 2:
            multi_copy += 1
        elif "Fragmented" in statuses:
            fragmented += 1
        # else: Missing -- no counter needed, derived below.

    complete = single_copy + multi_copy
    missing = n_markers - complete - fragmented

    def pct(n):
        return round((n / n_markers) * 100, 1) if n_markers else 0.0

    results = {
        "Complete percentage": pct(complete),
        "Single copy percentage": pct(single_copy),
        "Multi copy percentage": pct(multi_copy),
        "Fragmented percentage": pct(fragmented),
        "Missing percentage": pct(missing),
        "n_markers": n_markers,
        "Complete": complete,
        "Multi copy": multi_copy,
        "dataset": args.lineage,
    }

    import json
    with open(args.output_json, "w") as handle:
        json.dump({"results": results}, handle, indent=2)

    prefix = f"{args.sample_id}.{args.stage}.{args.lineage}"
    one_line = (
        f"C:{results['Complete percentage']}%"
        f"[S:{results['Single copy percentage']}%,D:{results['Multi copy percentage']}%],"
        f"F:{results['Fragmented percentage']}%,M:{results['Missing percentage']}%,n:{n_markers}\n"
    )
    with open(args.output_txt, "w") as handle:
        handle.write(f"# BUSCO version is: merged from {len(args.full_tables)} chunk(s)\n")
        handle.write(f"# The lineage dataset is: {args.lineage} (number of BUSCOs: {n_markers})\n")
        handle.write(f"# Summarized benchmarking in BUSCO notation for {prefix}\n")
        handle.write(f"# BUSCO was run in mode: genome (split into {len(args.full_tables)} chunk(s) and merged)\n")
        handle.write("\n")
        handle.write(f"\t{one_line}")
        handle.write(f"\t{complete}\tComplete BUSCOs (C)\n")
        handle.write(f"\t{single_copy}\tComplete and single-copy BUSCOs (S)\n")
        handle.write(f"\t{multi_copy}\tComplete and duplicated BUSCOs (D)\n")
        handle.write(f"\t{fragmented}\tFragmented BUSCOs (F)\n")
        handle.write(f"\t{missing}\tMissing BUSCOs (M)\n")
        handle.write(f"\t{n_markers}\tTotal BUSCO groups searched\n")


if __name__ == "__main__":
    main()
