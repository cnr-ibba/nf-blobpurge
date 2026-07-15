#!/usr/bin/env python3
"""Build a comparative BUSCO completeness/duplication table across stages and lineages.

Reads BUSCO's own `short_summary...json` for every (stage, lineage) run and
writes a single TSV + JSON table, plus a per-lineage duplication-drop summary
between the 'filtered' stage and each purge stage. All numbers come straight
from BUSCO's JSON output -- nothing here is hand-entered.

Each input file must follow the naming convention produced by the BUSCO
module: `<sample_id>.<stage>.<lineage>.short_summary.json`, where stage is
one of raw/filtered/purge_dups/purge_haplotigs.
"""
import argparse
import json
import re
import sys

STAGE_ORDER = ["raw", "filtered", "purge_dups", "purge_haplotigs"]
FILENAME_RE = re.compile(r"^(?P<sample>.+)\.(?P<stage>raw|filtered|purge_dups|purge_haplotigs)\.(?P<lineage>.+)\.short_summary\.json$")


def load_busco_json(path):
    with open(path) as handle:
        data = json.load(handle)
    results = data.get("results", {})

    def pct(key):
        val = results.get(key)
        return float(val) if val is not None else None

    return {
        "complete_pct": pct("Complete percentage"),
        "single_copy_pct": pct("Single copy percentage"),
        "duplicated_pct": pct("Multi copy percentage"),
        "fragmented_pct": pct("Fragmented percentage"),
        "missing_pct": pct("Missing percentage"),
        "n_markers": results.get("n_markers"),
        "complete_n": results.get("Complete"),
        "duplicated_n": results.get("Multi copy"),
        "dataset": results.get("dataset") or data.get("lineage_dataset", {}).get("name"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("summaries", nargs="+", help="BUSCO short_summary*.json files, named <sample>.<stage>.<lineage>.short_summary.json")
    parser.add_argument("--dup-drop-threshold", type=float, default=0.5)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    rows = []
    for path in args.summaries:
        match = FILENAME_RE.match(path.split("/")[-1])
        if not match:
            sys.exit(f"ERROR: '{path}' does not match the expected BUSCO summary naming convention <sample>.<stage>.<lineage>.short_summary.json")
        parsed = load_busco_json(path)
        rows.append({"sample_id": args.sample_id, "stage": match.group("stage"), "lineage": match.group("lineage"), **parsed})

    rows.sort(key=lambda r: (r["lineage"], STAGE_ORDER.index(r["stage"]) if r["stage"] in STAGE_ORDER else 99))

    fieldnames = [
        "sample_id", "stage", "lineage", "dataset", "n_markers",
        "complete_pct", "single_copy_pct", "duplicated_pct", "fragmented_pct", "missing_pct",
        "complete_n", "duplicated_n",
    ]
    with open(args.output_tsv, "w") as handle:
        handle.write("\t".join(fieldnames) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(f, "")) for f in fieldnames) + "\n")

    # Duplication drop: filtered -> each later purge stage, per lineage.
    by_lineage_stage = {}
    for row in rows:
        by_lineage_stage.setdefault(row["lineage"], {})[row["stage"]] = row

    duplication_drop = {}
    for lineage, stages in by_lineage_stage.items():
        filtered_dup = stages.get("filtered", {}).get("duplicated_pct")
        duplication_drop[lineage] = {}
        for purge_stage in ("purge_dups", "purge_haplotigs"):
            purge_dup = stages.get(purge_stage, {}).get("duplicated_pct")
            if filtered_dup is None or purge_dup is None:
                duplication_drop[lineage][purge_stage] = None
                continue
            if filtered_dup == 0:
                relative_drop = 0.0
            else:
                relative_drop = (filtered_dup - purge_dup) / filtered_dup
            duplication_drop[lineage][purge_stage] = {
                "filtered_duplicated_pct": filtered_dup,
                "purged_duplicated_pct": purge_dup,
                "relative_drop": relative_drop,
                "meets_threshold": relative_drop >= args.dup_drop_threshold,
            }

    output = {
        "sample_id": args.sample_id,
        "rows": rows,
        "duplication_drop": duplication_drop,
        "dup_drop_threshold": args.dup_drop_threshold,
    }
    with open(args.output_json, "w") as handle:
        json.dump(output, handle, indent=2)


if __name__ == "__main__":
    main()
