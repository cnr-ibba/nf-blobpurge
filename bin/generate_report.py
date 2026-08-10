#!/usr/bin/env python3
"""Generate the per-sample blobpurge HTML report.

Every number in this report is read from an upstream process's own output
file (assembly_stats.json per stage, the BTK_FILTER span-check JSON, BUSCO's
comparative JSON, purge_dups/purge_haplotigs outputs) -- nothing is
hand-entered. Missing optional inputs (GenomeScope2 summary) are reported as
explicitly "not available", never silently dropped.
"""
import argparse
import html
import json
import re
import sys

STATS_FILENAME_RE = re.compile(r"^.+\.(?P<stage>raw|filtered|purge_dups|purge_haplotigs)\.stats\.json$")


def load_json(path):
    with open(path) as handle:
        return json.load(handle)


def load_stats(paths):
    """paths: ASSEMBLY_STATS outputs named <sample>.<stage>.stats.json -> {stage: stats_dict}"""
    stats = {}
    for path in paths or []:
        match = STATS_FILENAME_RE.match(path.split("/")[-1])
        if not match:
            sys.exit(f"ERROR: '{path}' does not match the expected assembly-stats naming convention <sample>.<stage>.stats.json")
        stats[match.group("stage")] = load_json(path)
    return stats


def parse_genomescope_summary(path):
    """Extract the haploid genome length estimate from a GenomeScope2 summary.txt."""
    if path is None:
        return None
    with open(path) as handle:
        text = handle.read()
    match = re.search(r"Genome Haploid Length\s+([\d,]+)\s*bp\s+([\d,]+)\s*bp", text)
    if not match:
        return {"raw_text": text, "haploid_length_bp": None}
    low = int(match.group(1).replace(",", ""))
    high = int(match.group(2).replace(",", ""))
    return {"raw_text": text, "haploid_length_bp": (low + high) // 2, "haploid_length_range": [low, high]}


def bar_svg(rows, value_key, label_key, max_value, width=420, bar_height=22, color="#2f6fed"):
    """Minimal, self-contained inline SVG horizontal bar chart (no external deps)."""
    svg_rows = []
    y = 0
    for row in rows:
        value = row[value_key] if row[value_key] is not None else 0
        frac = (value / max_value) if max_value else 0
        bar_w = max(2, frac * (width - 160))
        svg_rows.append(
            f'<text x="0" y="{y + bar_height * 0.7:.0f}" font-size="12">{html.escape(str(row[label_key]))}</text>'
            f'<rect x="160" y="{y + 3}" width="{bar_w:.1f}" height="{bar_height - 6}" fill="{color}" rx="3"/>'
            f'<text x="{160 + bar_w + 6:.1f}" y="{y + bar_height * 0.7:.0f}" font-size="12">'
            f'{value:.2f}%</text>'
        )
        y += bar_height
    return f'<svg width="{width}" height="{y}" xmlns="http://www.w3.org/2000/svg" font-family="sans-serif" fill="#1a1a1a">{"".join(svg_rows)}</svg>'


def fmt_bp(n):
    if n is None:
        return "n/a"
    return f"{n:,} bp"


def fmt_pct(x, digits=1):
    if x is None:
        return "n/a"
    return f"{x * 100:.{digits}f}%"


def compute_verdict(stats, busco_comparison, disagreement, dup_drop_threshold):
    span_filtered = stats.get("filtered", {}).get("total_span")
    span_pd = stats.get("purge_dups", {}).get("total_span")
    span_ph = stats.get("purge_haplotigs", {}).get("total_span")

    span_shrink_pd = None
    if span_filtered and span_pd:
        span_shrink_pd = (span_filtered - span_pd) / span_filtered

    dup_drops = busco_comparison.get("duplication_drop", {})
    pd_drops = [v["purge_dups"]["relative_drop"] for v in dup_drops.values() if v.get("purge_dups")]
    avg_pd_drop = sum(pd_drops) / len(pd_drops) if pd_drops else None

    meets_dup_threshold = avg_pd_drop is not None and avg_pd_drop >= dup_drop_threshold
    meets_span_shrink = span_shrink_pd is not None and span_shrink_pd > 0
    disagree_flag = disagreement is not None and disagreement.get("flagged", False)

    if meets_dup_threshold and meets_span_shrink and not disagree_flag:
        verdict = "YES"
        explanation = (
            "The BUSCO duplication reduction and the assembly size contraction "
            "after purging are consistent with uncollapsed heterozygous haplotigs "
            "as the main cause of the surplus."
        )
    elif (avg_pd_drop and avg_pd_drop > 0) or meets_span_shrink or disagree_flag:
        verdict = "PARTIAL"
        reasons = []
        if not meets_dup_threshold:
            reasons.append(f"the BUSCO duplication drop ({fmt_pct(avg_pd_drop)}) does not reach the threshold ({fmt_pct(dup_drop_threshold)})")
        if not meets_span_shrink:
            reasons.append("the assembly size does not shrink markedly after purging")
        if disagree_flag:
            reasons.append("purge_dups and purge_haplotigs diverge substantially")
        explanation = "Only partially concordant signals: " + "; ".join(reasons) + "."
    else:
        verdict = "NO"
        explanation = (
            "Neither the BUSCO duplication drop nor the assembly size contraction "
            "support the hypothesis of uncollapsed heterozygous haplotigs as the "
            "main cause of the surplus; consider other sources (residual "
            "contamination, repeats, assembly heterogeneity)."
        )

    return {
        "verdict": verdict,
        "explanation": explanation,
        "avg_purge_dups_dup_drop": avg_pd_drop,
        "span_shrink_purge_dups": span_shrink_pd,
        "disagreement_flagged": disagree_flag,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stats", nargs="+", required=True, metavar="PATH", help="ASSEMBLY_STATS outputs, named <sample>.<stage>.stats.json")
    parser.add_argument("--span-check", required=True)
    parser.add_argument("--busco-comparison", required=True)
    parser.add_argument("--genomescope-summary", default=None)
    parser.add_argument("--organelle-report", default=None, help="bin/detect_organelles.py --output-json, if --run_organelle_isolation was used")
    parser.add_argument("--caveats", action="append", default=[])
    parser.add_argument("--run-purge-haplotigs", action="store_true")
    parser.add_argument("--dup-drop-threshold", type=float, default=0.5)
    parser.add_argument("--disagreement-threshold", type=float, default=0.15)
    parser.add_argument("--output-html", required=True)
    args = parser.parse_args()

    stats = load_stats(args.stats)
    span_check = load_json(args.span_check)
    busco_comparison = load_json(args.busco_comparison)
    genomescope = parse_genomescope_summary(args.genomescope_summary)
    organelle_report = load_json(args.organelle_report) if args.organelle_report else None

    span_pd = stats.get("purge_dups", {}).get("total_span")
    span_ph = stats.get("purge_haplotigs", {}).get("total_span")
    disagreement = None
    if span_pd and span_ph:
        mean_span = (span_pd + span_ph) / 2
        rel_diff = abs(span_pd - span_ph) / mean_span
        disagreement = {
            "purge_dups_span": span_pd,
            "purge_haplotigs_span": span_ph,
            "relative_diff": rel_diff,
            "flagged": rel_diff > args.disagreement_threshold,
        }

    verdict = compute_verdict(stats, busco_comparison, disagreement, args.dup_drop_threshold)

    stage_labels = {"raw": "Original", "filtered": "Filtered (BTK_FILTER)", "purge_dups": "Purged (purge_dups)", "purge_haplotigs": "Purged (purge_haplotigs)"}
    stage_order = [s for s in ("raw", "filtered", "purge_dups", "purge_haplotigs") if s in stats]

    span_rows = "".join(
        f"<tr><td>{stage_labels.get(s, s)}</td><td>{fmt_bp(stats[s]['total_span'])}</td>"
        f"<td>{stats[s]['n_contigs']:,}</td><td>{fmt_bp(stats[s].get('n50'))}</td></tr>"
        for s in stage_order
    )

    busco_rows_by_lineage = {}
    for row in busco_comparison.get("rows", []):
        busco_rows_by_lineage.setdefault(row["lineage"], []).append(row)

    busco_sections = []
    for lineage, rows in busco_rows_by_lineage.items():
        chart_rows = [
            {"label": stage_labels.get(r["stage"], r["stage"]), "value": r["duplicated_pct"]}
            for r in sorted(rows, key=lambda r: stage_order.index(r["stage"]) if r["stage"] in stage_order else 99)
        ]
        table_rows = "".join(
            f"<tr><td>{stage_labels.get(r['stage'], r['stage'])}</td>"
            f"<td>{r['complete_pct']:.1f}%</td><td>{r['single_copy_pct']:.1f}%</td>"
            f"<td>{r['duplicated_pct']:.1f}%</td><td>{r['fragmented_pct']:.1f}%</td>"
            f"<td>{r['missing_pct']:.1f}%</td><td>{r['n_markers']}</td></tr>"
            for r in sorted(rows, key=lambda r: stage_order.index(r["stage"]) if r["stage"] in stage_order else 99)
        )
        busco_sections.append(f"""
        <h3>Lineage: {html.escape(lineage)}</h3>
        <div class="chart">{bar_svg(chart_rows, 'value', 'label', 100, color='#d1495b')}</div>
        <table>
        <thead><tr><th>Stage</th><th>Complete</th><th>Single</th><th>Duplicated</th><th>Fragmented</th><th>Missing</th><th># markers</th></tr></thead>
        <tbody>{table_rows}</tbody>
        </table>
        """)

    if genomescope is None:
        genomescope_html = '<p class="warn">GenomeScope2 not provided (--genomescope_summary absent): comparison with the k-mer estimate is NOT available.</p>'
    else:
        hap_len = genomescope.get("haploid_length_bp")
        if hap_len is None:
            genomescope_html = '<p class="warn">GenomeScope2 provided, but the "Genome Haploid Length" estimate was not found in summary.txt.</p>'
        else:
            final_span = stats.get("purge_dups", stats.get("filtered", {})).get("total_span")
            ratio = final_span / hap_len if final_span and hap_len else None
            genomescope_html = (
                f"<p>GenomeScope2 estimate (haploid): <b>{fmt_bp(hap_len)}</b> "
                f"(range {fmt_bp(genomescope['haploid_length_range'][0])} - {fmt_bp(genomescope['haploid_length_range'][1])}).<br/>"
                f"Purged assembly (purge_dups): <b>{fmt_bp(final_span)}</b> "
                f"({'+' if ratio and ratio > 1 else ''}{fmt_pct((ratio - 1)) if ratio else 'n/a'} relative to the estimate).</p>"
            )

    disagreement_html = '<p class="warn">purge_haplotigs not run: no cross-check available.</p>'
    if disagreement is not None:
        flag = disagreement["flagged"]
        css_class = "warn" if flag else "ok"
        disagreement_html = (
            f'<p class="{css_class}">purge_dups: {fmt_bp(disagreement["purge_dups_span"])} vs '
            f'purge_haplotigs: {fmt_bp(disagreement["purge_haplotigs_span"])} '
            f'(relative difference {fmt_pct(disagreement["relative_diff"])}, threshold {fmt_pct(args.disagreement_threshold)}). '
        )
        if flag:
            disagreement_html += (
                "The two tools diverge substantially: this should be resolved with "
                "long-read data, not arbitrated automatically by this pipeline."
            )
        else:
            disagreement_html += "The two tools agree within the configured threshold."
        disagreement_html += "</p>"

    caveats_html = "".join(f"<li>{html.escape(c)}</li>" for c in args.caveats) or "<li>No caveats.</li>"

    span_check_class = "ok" if span_check["pass"] else "fail"
    span_check_html = (
        f'<p class="{span_check_class}">span(filtered) + span(excluded from BlobDir) = '
        f'{fmt_bp(span_check["filtered_span"])} + {fmt_bp(span_check["excluded_span_from_blobdir"])} = '
        f'{fmt_bp(span_check["reconstructed_span"])}, vs span(original) {fmt_bp(span_check["original_span"])} '
        f'(relative diff {fmt_pct(span_check["relative_diff"])}, tolerance {fmt_pct(span_check["tolerance"])}). '
        f'{"OK" if span_check["pass"] else "FAILED"}</p>'
    )

    if organelle_report is None:
        organelle_html = (
            '<p class="warn">Organelle isolation not enabled (see --run_organelle_isolation): no separate '
            "accounting for mitochondrial/plastid contigs is available.</p>"
        )
    else:
        n_isolated = organelle_report.get("n_isolated", 0)
        threshold_note = (
            f'coverage &ge; {organelle_report.get("threshold", 0):.1f}x '
            f'({organelle_report.get("coverage_multiplier")}x the '
            f'{organelle_report.get("baseline_coverage", 0):.1f}x nuclear baseline)'
        )
        if n_isolated:
            reference_used = organelle_report.get("reference_fasta_used", False)
            ref_col = "<th>Reference match</th>" if reference_used else ""
            contig_rows = "".join(
                f"<tr><td>{html.escape(c['id'])}</td><td>{fmt_bp(c['length'])}</td>"
                f"<td>{c['mean_depth']:.1f}x</td><td>{(c['coverage_ratio'] or 0):.1f}x</td>"
                + (f"<td>{'Yes' if c.get('reference_matched') else 'No'}</td>" if reference_used else "")
                + "</tr>"
                for c in organelle_report.get("contigs", [])
                if c.get("isolated")
            )
            organelle_html = (
                f'<p class="ok">{n_isolated} organelle-like contig(s) isolated '
                f'(total {fmt_bp(organelle_report.get("total_isolated_span"))}), {threshold_note}. '
                "Excluded from purge_dups/purge_haplotigs coverage-cutoff estimation and re-merged into the "
                "final assembly at each stage. A matching paired-end read subset was exported alongside the "
                "isolated FASTA (see the pipeline output directory) for use with external organelle-assembly "
                "tools such as GetOrganelle, MitoHiFi, or oatk.</p>"
                f"<table><thead><tr><th>Contig</th><th>Length</th><th>Mean depth</th><th>Coverage ratio</th>{ref_col}</tr></thead>"
                f"<tbody>{contig_rows}</tbody></table>"
            )
        else:
            organelle_html = f'<p class="warn">Organelle isolation enabled, but no contig was isolated ({threshold_note}).</p>'

    verdict_class = {"YES": "ok", "PARTIAL": "warn", "NO": "fail"}[verdict["verdict"]]

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>blobpurge report - {html.escape(args.sample_id)}</title>
<style>
  body {{ font-family: -apple-system, Helvetica, Arial, sans-serif; max-width: 960px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; background: #fff; }}
  h1 {{ border-bottom: 3px solid #2f6fed; padding-bottom: .5rem; }}
  h2 {{ margin-top: 2.5rem; border-bottom: 1px solid #ddd; padding-bottom: .3rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: .75rem 0 1.5rem; font-size: .92rem; }}
  th, td {{ border: 1px solid #ddd; padding: .4rem .6rem; text-align: right; }}
  th:first-child, td:first-child {{ text-align: left; }}
  th {{ background: #f4f6fb; }}
  .ok {{ color: #1a7f37; }}
  .warn {{ color: #9a6700; }}
  .fail {{ color: #cf222e; font-weight: bold; }}
  .verdict-box {{ border: 2px solid currentColor; border-radius: 8px; padding: 1rem 1.25rem; font-size: 1.1rem; margin: 1rem 0 2rem; }}
  .caveats {{ background: #fff8e6; border-left: 4px solid #9a6700; padding: .75rem 1rem; }}
  .chart {{ margin: .5rem 0; }}
  code {{ background: #f4f6fb; padding: .1rem .3rem; border-radius: 3px; }}
</style>
</head>
<body>
<h1>blobpurge report &mdash; {html.escape(args.sample_id)}</h1>

<h2>Verdict</h2>
<div class="verdict-box {verdict_class}">
  <b>{verdict['verdict']}</b> &mdash; {html.escape(verdict['explanation'])}
</div>

<h2>Assembly span per stage</h2>
<table>
<thead><tr><th>Stage</th><th>Total span</th><th># contigs</th><th>N50</th></tr></thead>
<tbody>{span_rows}</tbody>
</table>

<h3>Span-conservation check (BTK_FILTER)</h3>
{span_check_html}

<h2>Organelle contig isolation</h2>
{organelle_html}

<h2>Comparison with GenomeScope2</h2>
{genomescope_html}

<h2>Comparative BUSCO</h2>
{''.join(busco_sections)}

<h2>purge_dups vs purge_haplotigs</h2>
{disagreement_html}

<h2>Notes and caveats</h2>
<div class="caveats"><ul>{caveats_html}</ul></div>

<p style="color:#666; font-size:.85rem; margin-top:3rem;">Automatically generated by the blobpurge pipeline. Every number above is read from upstream processes' own outputs (nothing is hand-entered).</p>
</body>
</html>
"""

    with open(args.output_html, "w") as handle:
        handle.write(html_doc)

    summary = {
        "sample_id": args.sample_id,
        "verdict": verdict,
        "span_check": span_check,
        "disagreement": disagreement,
        "genomescope_provided": genomescope is not None,
        "organelle": organelle_report,
    }
    with open(args.output_html.replace(".html", ".summary.json"), "w") as handle:
        json.dump(summary, handle, indent=2)


if __name__ == "__main__":
    main()
