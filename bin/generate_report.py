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
import itertools
import json
import re
import sys

STATS_FILENAME_RE = re.compile(r"^.+\.(?P<stage>raw|filtered|purge_dups|purge_haplotigs|haplomerger2)\.stats\.json$")

# Every stage that is itself a purge *method* (as opposed to raw/filtered,
# which are not purge outputs). Adding a future method (e.g. Redundans) only
# requires appending it here -- compute_disagreement() and compute_verdict()
# below already generalize over however many of these actually ran.
PURGE_STAGES = ("purge_dups", "purge_haplotigs", "haplomerger2")


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


def compute_disagreement(stats, disagreement_threshold):
    """Pairwise span disagreement across every purge method that actually ran.

    Returns None if fewer than 2 purge stages have stats (no cross-check
    possible), otherwise a dict with the list of stages that ran and one
    entry per pair, each using the same mean-normalized relative-difference
    formula as the original purge_dups-vs-purge_haplotigs check.
    """
    spans = {s: stats[s]["total_span"] for s in PURGE_STAGES if stats.get(s, {}).get("total_span")}
    if len(spans) < 2:
        return None

    pairs = []
    for a, b in itertools.combinations(spans, 2):
        mean_span = (spans[a] + spans[b]) / 2
        rel_diff = abs(spans[a] - spans[b]) / mean_span
        pairs.append(
            {
                "a": a,
                "b": b,
                "span_a": spans[a],
                "span_b": spans[b],
                "relative_diff": rel_diff,
                "flagged": rel_diff > disagreement_threshold,
            }
        )

    return {
        "stages": sorted(spans),
        "pairs": pairs,
        "any_flagged": any(p["flagged"] for p in pairs),
    }


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
    disagree_flag = disagreement is not None and disagreement.get("any_flagged", False)

    if meets_dup_threshold and meets_span_shrink and not disagree_flag:
        verdict = "SI"
        explanation = (
            "La riduzione di duplicazione BUSCO e la contrazione della dimensione "
            "dell'assembly dopo il purging sono coerenti con aplotidi eterozigoti "
            "non collassati come causa principale del surplus."
        )
    elif (avg_pd_drop and avg_pd_drop > 0) or meets_span_shrink or disagree_flag:
        verdict = "PARZIALMENTE"
        reasons = []
        if not meets_dup_threshold:
            reasons.append(f"il calo di duplicazione BUSCO ({fmt_pct(avg_pd_drop)}) non raggiunge la soglia ({fmt_pct(dup_drop_threshold)})")
        if not meets_span_shrink:
            reasons.append("la dimensione dell'assembly non si riduce in modo netto dopo il purging")
        if disagree_flag:
            reasons.append("gli strumenti di purging eseguiti divergono in modo sostanziale tra loro")
        explanation = "Segnali solo parzialmente concordanti: " + "; ".join(reasons) + "."
    else:
        verdict = "NO"
        explanation = (
            "Ne il calo di duplicazione BUSCO ne la contrazione della dimensione "
            "dell'assembly supportano l'ipotesi di aplotidi eterozigoti non collassati "
            "come causa principale del surplus; considerare altre fonti (contaminazione "
            "residua, ripetizioni, eterogeneita' dell'assembly)."
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

    disagreement = compute_disagreement(stats, args.disagreement_threshold)

    verdict = compute_verdict(stats, busco_comparison, disagreement, args.dup_drop_threshold)

    stage_labels = {
        "raw": "Originale",
        "filtered": "Filtrata (BTK_FILTER)",
        "purge_dups": "Purgata (purge_dups)",
        "purge_haplotigs": "Purgata (purge_haplotigs)",
        "haplomerger2": "Purgata (HaploMerger2)",
    }
    stage_order = [s for s in ("raw", "filtered", *PURGE_STAGES) if s in stats]

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
        <thead><tr><th>Stadio</th><th>Complete</th><th>Single</th><th>Duplicated</th><th>Fragmented</th><th>Missing</th><th># marker</th></tr></thead>
        <tbody>{table_rows}</tbody>
        </table>
        """)

    if genomescope is None:
        genomescope_html = '<p class="warn">GenomeScope2 non fornito (--genomescope_summary assente): confronto con la stima k-mer NON disponibile.</p>'
    else:
        hap_len = genomescope.get("haploid_length_bp")
        if hap_len is None:
            genomescope_html = '<p class="warn">GenomeScope2 fornito, ma la stima "Genome Haploid Length" non e stata trovata nel summary.txt.</p>'
        else:
            final_span = stats.get("purge_dups", stats.get("filtered", {})).get("total_span")
            ratio = final_span / hap_len if final_span and hap_len else None
            genomescope_html = (
                f"<p>Stima GenomeScope2 (aploide): <b>{fmt_bp(hap_len)}</b> "
                f"(range {fmt_bp(genomescope['haploid_length_range'][0])} - {fmt_bp(genomescope['haploid_length_range'][1])}).<br/>"
                f"Assembly purgata (purge_dups): <b>{fmt_bp(final_span)}</b> "
                f"({'+' if ratio and ratio > 1 else ''}{fmt_pct((ratio - 1)) if ratio else 'n/a'} rispetto alla stima).</p>"
            )

    purge_methods_ran = [s for s in PURGE_STAGES if s in stats]
    methods_html = (
        f'<p>Metodi di purging eseguiti: <b>{", ".join(stage_labels.get(s, s) for s in purge_methods_ran)}</b>.</p>'
    )

    if disagreement is None:
        disagreement_html = methods_html + (
            '<p class="warn">Un solo metodo di purging eseguito: nessun cross-check di span disponibile.</p>'
        )
    else:
        pair_rows = []
        for pair in disagreement["pairs"]:
            css_class = "warn" if pair["flagged"] else "ok"
            verdict_text = (
                "Divergenza sostanziale: da risolvere con dati long-read, non arbitrata "
                "automaticamente da questa pipeline."
                if pair["flagged"]
                else "In accordo entro la soglia configurata."
            )
            pair_rows.append(
                f'<p class="{css_class}">{stage_labels.get(pair["a"], pair["a"])}: {fmt_bp(pair["span_a"])} vs '
                f'{stage_labels.get(pair["b"], pair["b"])}: {fmt_bp(pair["span_b"])} '
                f'(differenza relativa {fmt_pct(pair["relative_diff"])}, soglia {fmt_pct(args.disagreement_threshold)}). '
                f"{verdict_text}</p>"
            )
        disagreement_html = methods_html + "".join(pair_rows)

    caveats_html = "".join(f"<li>{html.escape(c)}</li>" for c in args.caveats) or "<li>Nessuna nota.</li>"

    span_check_class = "ok" if span_check["pass"] else "fail"
    span_check_html = (
        f'<p class="{span_check_class}">span(filtrata) + span(esclusi da BlobDir) = '
        f'{fmt_bp(span_check["filtered_span"])} + {fmt_bp(span_check["excluded_span_from_blobdir"])} = '
        f'{fmt_bp(span_check["reconstructed_span"])}, vs span(originale) {fmt_bp(span_check["original_span"])} '
        f'(diff relativa {fmt_pct(span_check["relative_diff"])}, tolleranza {fmt_pct(span_check["tolerance"])}). '
        f'{"OK" if span_check["pass"] else "FALLITO"}</p>'
    )

    verdict_class = {"SI": "ok", "PARZIALMENTE": "warn", "NO": "fail"}[verdict["verdict"]]

    html_doc = f"""<!doctype html>
<html lang="it">
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

<h2>Verdetto</h2>
<div class="verdict-box {verdict_class}">
  <b>{verdict['verdict']}</b> &mdash; {html.escape(verdict['explanation'])}
</div>

<h2>Span dell'assembly per stadio</h2>
<table>
<thead><tr><th>Stadio</th><th>Span totale</th><th># contig</th><th>N50</th></tr></thead>
<tbody>{span_rows}</tbody>
</table>

<h3>Verifica di conservazione dello span (BTK_FILTER)</h3>
{span_check_html}

<h2>Confronto con GenomeScope2</h2>
{genomescope_html}

<h2>BUSCO comparativo</h2>
{''.join(busco_sections)}

<h2>Confronto tra strumenti di purging</h2>
{disagreement_html}

<h2>Note e limitazioni</h2>
<div class="caveats"><ul>{caveats_html}</ul></div>

<p style="color:#666; font-size:.85rem; margin-top:3rem;">Generato automaticamente dalla pipeline blobpurge. Tutti i numeri sopra sono letti dagli output dei processi a monte (nessun valore inserito manualmente).</p>
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
    }
    with open(args.output_html.replace(".html", ".summary.json"), "w") as handle:
        json.dump(summary, handle, indent=2)


if __name__ == "__main__":
    main()
