#!/usr/bin/env python3
"""Split a FASTA into span-balanced chunks so BUSCO can run on each separately.

BUSCO's own single-threaded post-processing of miniprot/metaeuk candidate
alignments does not scale down with contig fragmentation -- it scales with
total genome content searched against the reference protein set, so it can
run out of memory on assemblies well within normal size ranges regardless of
how clean/unfragmented they are. Splitting the assembly into several
span-balanced chunks and running BUSCO on each independently keeps each
individual BUSCO invocation's memory bounded; results are merged back
per-ortholog downstream (see bin/merge_busco_chunks.py).

Contigs are never split mid-sequence: chunk assignment walks contigs in
FASTA input order and closes a chunk once adding the next contig would push
it over --max-length (a single contig longer than --max-length still becomes
its own oversized chunk). --max-length <= 0 means "unlimited" -> always
exactly one chunk.
"""
import argparse
import gzip
import json


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_lengths(fasta_path):
    """Pass 1: contig ids and lengths, in FASTA order. Sequence is not buffered."""
    lengths = []
    current_id = None
    current_len = 0
    with open_maybe_gzip(fasta_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths.append((current_id, current_len))
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id is not None:
            lengths.append((current_id, current_len))
    return lengths


def assign_chunks(lengths, max_span):
    """Sequential accumulation in input order. Returns {contig_id: chunk_index (1-based)}."""
    assignment = {}
    chunk_index = 1
    chunk_span = 0
    chunk_has_contig = False
    for contig_id, length in lengths:
        if max_span > 0 and chunk_has_contig and chunk_span + length > max_span:
            chunk_index += 1
            chunk_span = 0
            chunk_has_contig = False
        assignment[contig_id] = chunk_index
        chunk_span += length
        chunk_has_contig = True
    return assignment


def write_chunks(fasta_path, assignment, n_chunks, output_prefix):
    """Pass 2: stream the FASTA again, writing each record to its assigned chunk file."""
    handles = {
        k: open(f"{output_prefix}.chunk{k}of{n_chunks}.fasta", "w")
        for k in range(1, n_chunks + 1)
    }
    try:
        current_id = None
        current_handle = None
        with open_maybe_gzip(fasta_path) as handle:
            for line in handle:
                if line.startswith(">"):
                    current_id = line[1:].split()[0]
                    current_handle = handles[assignment[current_id]]
                current_handle.write(line if line.endswith("\n") else line + "\n")
    finally:
        for h in handles.values():
            h.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--max-span", type=int, required=True, help="Target max span (bp) per chunk; <=0 disables splitting")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    lengths = read_lengths(args.fasta)
    assignment = assign_chunks(lengths, args.max_span)
    n_chunks = max(assignment.values()) if assignment else 1

    if not lengths:
        # Degenerate but keep the pipeline moving with a single empty chunk
        # rather than special-casing an empty assembly downstream.
        open(f"{args.output_prefix}.chunk1of1.fasta", "w").close()
        n_chunks = 1

    if lengths:
        write_chunks(args.fasta, assignment, n_chunks, args.output_prefix)

    chunk_stats = {k: {"n_contigs": 0, "span": 0} for k in range(1, n_chunks + 1)}
    for contig_id, length in lengths:
        k = assignment[contig_id]
        chunk_stats[k]["n_contigs"] += 1
        chunk_stats[k]["span"] += length

    report = {
        "sample_id": args.sample_id,
        "stage": args.stage,
        "chunk_max_span": args.max_span,
        "n_chunks": n_chunks,
        "total_span": sum(length for _cid, length in lengths),
        "total_contigs": len(lengths),
        "chunks": [
            {"index": k, "n_contigs": chunk_stats[k]["n_contigs"], "span": chunk_stats[k]["span"]}
            for k in range(1, n_chunks + 1)
        ],
    }
    with open(args.output_json, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
