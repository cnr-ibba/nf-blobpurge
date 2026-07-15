#!/usr/bin/env python3
"""Deduce the taxrule used to build a BlobDir, from its meta.json.

blobtools filter --summary defaults to --taxrule bestsumorder, which does not
necessarily exist in every BlobDir (it depends on the taxrule(s) the BlobDir
was built with). This script inspects meta.json and prints a single taxrule
name to stdout, or exits non-zero with an explicit, actionable error.
"""
import argparse
import json
import sys

KNOWN_TAXRULE_PREFIXES = ("bestsumorder", "bestsum")


def candidate_taxrules_from_fields(meta):
    fields = meta.get("fields", [])
    found = set()
    for field in fields:
        field_id = field.get("id", "") if isinstance(field, dict) else ""
        for prefix in KNOWN_TAXRULE_PREFIXES:
            if field_id == f"{prefix}_phylum" or field_id.startswith(f"{prefix}_"):
                found.add(prefix)
    return sorted(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("meta_json", help="Path to <blobdir>/meta.json")
    args = parser.parse_args()

    with open(args.meta_json) as handle:
        meta = json.load(handle)

    # 1) Explicit settings in meta.json (most reliable, when present).
    settings = meta.get("settings", {})
    if isinstance(settings.get("taxrule"), str) and settings["taxrule"]:
        print(settings["taxrule"])
        return
    taxrules = settings.get("taxrules")
    if isinstance(taxrules, list) and taxrules:
        print(taxrules[0])
        return

    # 2) Infer from the category fields actually present in the BlobDir.
    candidates = candidate_taxrules_from_fields(meta)
    if len(candidates) == 1:
        print(candidates[0])
        return

    field_ids = sorted(
        f.get("id", "?") for f in meta.get("fields", []) if isinstance(f, dict)
    )
    if len(candidates) > 1:
        sys.exit(
            "ERROR: could not unambiguously deduce --taxrule from meta.json: "
            f"multiple candidate taxrules found ({', '.join(candidates)}). "
            "Pass --taxrule explicitly. Fields available in this BlobDir: "
            f"{', '.join(field_ids)}"
        )

    sys.exit(
        "ERROR: could not deduce --taxrule from this BlobDir's meta.json "
        "(no 'settings.taxrule(s)' entry and no bestsum/bestsumorder fields "
        "found). Pass --taxrule explicitly, matching one of the taxonomy "
        f"fields present in this BlobDir: {', '.join(field_ids)}"
    )


if __name__ == "__main__":
    main()
