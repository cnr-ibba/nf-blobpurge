#!/usr/bin/env python3
"""Drop contigs shorter than a minimum length from a FASTA file, as JSON.

Used to keep BUSCO's raw-stage input from ballooning on highly fragmented
assemblies: contigs too short to ever hold a complete gene model still get
searched by miniprot, and BUSCO's own single-threaded post-processing of
that candidate-alignment output can run out of memory as a result. Streams
one record at a time (never buffers the whole assembly) so this stays cheap
regardless of assembly size.
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
    parser.add_argument("--stage", required=True)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-mqc-json", required=True, help="MultiQC custom-content JSON path")
    args = parser.parse_args()

    n_contigs_in = 0
    n_contigs_out = 0
    span_in = 0
    span_out = 0

    header = None
    seq_lines = []
    seq_len = 0

    def flush(out_handle):
        nonlocal n_contigs_out, span_out
        if header is None:
            return
        if seq_len >= args.min_length:
            n_contigs_out += 1
            span_out += seq_len
            out_handle.write(header)
            out_handle.writelines(seq_lines)

    with open_maybe_gzip(args.fasta) as in_handle, open(args.output_fasta, "w") as out_handle:
        for line in in_handle:
            if line.startswith(">"):
                flush(out_handle)
                n_contigs_in += 1
                header = line
                seq_lines = []
                seq_len = 0
            else:
                seq_lines.append(line)
                seq_len += len(line.strip())
                span_in += len(line.strip())
        flush(out_handle)

    n_contigs_removed = n_contigs_in - n_contigs_out
    bp_removed = span_in - span_out

    report = {
        "sample_id": args.sample_id,
        "stage": args.stage,
        "min_contig_length": args.min_length,
        "n_contigs_in": n_contigs_in,
        "n_contigs_out": n_contigs_out,
        "n_contigs_removed": n_contigs_removed,
        "span_in": span_in,
        "span_out": span_out,
        "bp_removed": bp_removed,
        "pct_contigs_removed": (100 * n_contigs_removed / n_contigs_in) if n_contigs_in else 0.0,
        "pct_bp_removed": (100 * bp_removed / span_in) if span_in else 0.0,
    }
    with open(args.output_json, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))

    mqc = {
        "id": "contig_length_filter",
        "section_name": "BUSCO raw-stage contig-length filter",
        "description": (
            "Contigs shorter than --busco_min_contig_length excluded from BUSCO's "
            "raw-stage input only (assembly span/N50 elsewhere are unaffected)."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "contig_length_filter_table",
            "title": "BUSCO raw-stage contig-length filter",
            "namespace": "Contig Length Filter",
        },
        "data": {
            f"{args.sample_id}_{args.stage}": {
                "Sample": args.sample_id,
                "Stage": args.stage,
                "Threshold (bp)": args.min_length,
                "Contigs in": report["n_contigs_in"],
                "Contigs removed": report["n_contigs_removed"],
                "bp removed": report["bp_removed"],
                "% bp removed": round(report["pct_bp_removed"], 2),
            }
        },
    }
    with open(args.output_mqc_json, "w") as handle:
        json.dump(mqc, handle, indent=2)


if __name__ == "__main__":
    main()
