#!/usr/bin/env python3
"""Generate the tiny synthetic dataset used by -profile test.

Produces:
  - assembly.fasta: 7 contigs -- a heterozygous duplicate pair (ctg1/ctg2,
    8% divergent, so BTK_FILTER->purge_dups has something real to purge),
    two unrelated unique host contigs, two "contaminant" contigs, and one
    organelle-like contig at extreme coverage.
  - blobdir/: a minimal, hand-built BlobDir (meta.json + identifiers/length/
    buscoregions_phylum field JSON) matching the schema BTK_FILTER expects.
  - reads.cram(.crai): paired-end alignments sampled independently from ctg1
    and ctg2 (each at half the depth used for the unique host contigs) --
    reflecting the real biology of an uncollapsed heterozygous locus, where
    each assembled copy only receives reads from the haplotype it truly
    represents. The pipeline only accepts a pre-aligned reads CRAM (the
    sanger-tol/blobtoolkit precondition -- see read_coverage.nf), so this
    fixture hand-builds the alignment directly instead of simulating FASTQ
    and running a real aligner: each pair's true contig-of-origin and
    fragment position are already known from sampling, so the ground-truth
    SAM record is a trivial full-length match (no aligner-introduced
    ambiguity to reason about, unlike the old bwa-mem2-based design this
    replaced). Sequencing errors are still injected via mutate() to keep
    per-base identity non-trivial for downstream coverage/GC stats. The
    organelle-like contig is sampled at 8x the normal host depth (its real
    biological signature -- high copy number per cell). No reads are
    simulated for the contaminant contigs.
  - organelle_reference.fasta: a 2%-divergent copy of the organelle contig,
    for --organelle_reference_fasta's minimap2 corroboration path.
  - samplesheet.csv

The organelle-like contig is deliberately labelled with the same BlobDir
taxon ("Pseudomonadota") as one of the "contaminant" contigs -- already in
conf/test.config's --exclude_taxa -- rather than a legitimate host taxon.
This is what actually exercises ORGANELLE_ISOLATE running *before*
BTK_FILTER: with --run_organelle_isolation true, the organelle contig must
survive into every downstream stage FASTA despite matching an excluded
taxon (were it labelled as host taxon instead, BTK_FILTER would never have
touched it anyway, and the fixture would not catch a regression back to the
old, buggy post-filter ordering).

Everything here is synthetic and deterministic (fixed RNG seed) -- this is a
plumbing fixture for -profile test, not a biologically meaningful assembly.
Building the CRAM shells out to `samtools` via Docker (not required to be on
PATH), using the same image CRAM_TO_BAM runs in.
"""
import gzip
import json
import os
import random
import subprocess
import tempfile
from pathlib import Path

random.seed(42)
OUT = Path(__file__).parent
BASES = "ACGT"
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.24--h9dcdb79_0"


def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def write_blobdir_json(path, data, *, gzip_compress=False):
    """Write a BlobDir field JSON file, optionally gzip-compressed
    (<path>.gz) to mirror newer BlobToolKit output. Compression is
    deliberately mixed across this fixture's files so -profile test
    exercises both the "prefer .json.gz" and "fall back to plain .json"
    branches of bin/detect_taxrule.py and bin/check_span.py's BlobDir-file
    resolution in a single run.
    """
    text = json.dumps(data, indent=2)
    if gzip_compress:
        with gzip.open(f"{path}.gz", "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def mutate(seq, rate):
    bases = list(seq)
    for i in range(len(bases)):
        if random.random() < rate:
            bases[i] = random.choice([b for b in BASES if b != bases[i]])
    return "".join(bases)


def base_composition(seq):
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    n = len(seq) - at - gc
    return round(gc / (at + gc), 4), n


def revcomp(seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def write_fasta(path, records):
    with open(path, "w") as handle:
        for rid, seq in records.items():
            handle.write(f">{rid}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")


def sample_reads(seq, ref_id, name_prefix, depth, read_len=100, insert=300, error_rate=0.005):
    """Sample paired-end fragments and return their ground-truth alignment:
    each pair already knows its true reference contig and 0-based leftmost
    mapping position for R1/R2, since it was sampled directly from `seq`
    (only substitution errors are injected -- no indels -- so every read's
    true CIGAR is a trivial full-length match). Both SEQ fields are kept in
    reference/forward-strand orientation (SAM convention for a read mapped
    to the reverse strand), not sequencer-readout orientation.
    """
    n_pairs = max(1, (len(seq) * depth) // (2 * read_len))
    pairs = []
    for i in range(n_pairs):
        if len(seq) <= insert:
            start = 0
        else:
            start = random.randint(0, len(seq) - insert)
        frag = seq[start:start + insert]
        if len(frag) < read_len:
            frag = (frag + seq)[:read_len * 2]
        pairs.append({
            "name": f"{name_prefix}_{i}",
            "ref": ref_id,
            "r1_pos": start,
            "r1_seq": mutate(frag[:read_len], error_rate),
            "r2_pos": start + len(frag) - read_len,
            "r2_seq": mutate(frag[-read_len:], error_rate),
        })
    return pairs


def write_cram(path, pairs, ref_lengths, read_len=100):
    """Hand-build a SAM (header + one properly-paired record pair per
    fragment, MAPQ 60, full-length M CIGAR -- see sample_reads) and convert
    it to a coordinate-sorted, indexed CRAM against ref_lengths' contig
    order/lengths via a one-off samtools container. This mirrors what
    sanger-tol/blobtoolkit's own upstream alignment step produces: a reads
    CRAM aligned to the *raw*, unfiltered assembly.
    """
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    for ref_id, length in ref_lengths.items():
        lines.append(f"@SQ\tSN:{ref_id}\tLN:{length}")

    qual = "I" * read_len
    for pair in pairs:
        r1_end = pair["r1_pos"] + read_len
        r2_end = pair["r2_pos"] + read_len
        tlen = max(r1_end, r2_end) - min(pair["r1_pos"], pair["r2_pos"])
        lines.append("\t".join([
            pair["name"], "99", pair["ref"], str(pair["r1_pos"] + 1), "60",
            f"{read_len}M", "=", str(pair["r2_pos"] + 1), str(tlen),
            pair["r1_seq"], qual,
        ]))
        lines.append("\t".join([
            pair["name"], "147", pair["ref"], str(pair["r2_pos"] + 1), "60",
            f"{read_len}M", "=", str(pair["r1_pos"] + 1), str(-tlen),
            pair["r2_seq"], qual,
        ]))

    with tempfile.NamedTemporaryFile("w", suffix=".sam", dir=OUT, delete=False) as handle:
        handle.write("\n".join(lines) + "\n")
        sam_path = Path(handle.name)
    try:
        run_samtools([
            "sh", "-c",
            f"samtools sort -O cram --reference assembly.fasta "
            f"-o {path.name} {sam_path.name} && samtools index {path.name}",
        ])
    finally:
        sam_path.unlink(missing_ok=True)


def run_samtools(cmd):
    subprocess.run(
        [
            "docker", "run", "--rm",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "-v", f"{OUT}:/data", "-w", "/data",
            SAMTOOLS_IMAGE,
            *cmd,
        ],
        check=True,
    )


def main():
    ctg1 = random_seq(3000)
    ctg2 = mutate(ctg1, 0.08)
    # Unique host contigs are kept substantially larger than the duplicated
    # locus so the assembly-wide depth histogram's tallest peak is the
    # normal single-copy coverage level, not the half-depth haplotig region
    # -- matching what purge_haplotigs' cutoff estimation assumes (most of a
    # real assembly is single-copy; only a minority is uncollapsed haplotig).
    ctg3 = random_seq(12000)
    ctg6 = random_seq(8000)
    ctg4 = random_seq(1200)
    ctg5 = random_seq(1000)
    ctg7 = random_seq(1500)

    records = {
        "ctg1_host_A": ctg1,
        "ctg2_host_A_hap": ctg2,
        "ctg3_host_B": ctg3,
        "ctg4_contam_bac1": ctg4,
        "ctg5_contam_bac2": ctg5,
        "ctg6_host_C": ctg6,
        "ctg7_organelle": ctg7,
    }
    write_fasta(OUT / "assembly.fasta", records)
    write_fasta(OUT / "organelle_reference.fasta", {"organelle_ref": mutate(ctg7, 0.02)})

    depth = 80
    all_pairs = []
    all_pairs += sample_reads(ctg1, "ctg1_host_A", "h1", depth=depth // 2)
    all_pairs += sample_reads(ctg2, "ctg2_host_A_hap", "h2", depth=depth // 2)
    all_pairs += sample_reads(ctg3, "ctg3_host_B", "ctg3", depth=depth)
    all_pairs += sample_reads(ctg6, "ctg6_host_C", "ctg6", depth=depth)
    all_pairs += sample_reads(ctg7, "ctg7_organelle", "ctg7", depth=depth * 8)

    write_cram(OUT / "reads.cram", all_pairs, {rid: len(seq) for rid, seq in records.items()})

    # --- BlobDir ---
    blobdir = OUT / "blobdir"
    blobdir.mkdir(exist_ok=True)

    identifiers = list(records.keys())
    lengths = [len(records[i]) for i in identifiers]
    gc_ncount = [base_composition(records[i]) for i in identifiers]
    gc_portions = [gc for gc, _n in gc_ncount]
    n_counts = [n for _gc, n in gc_ncount]
    keys = ["Chlorophyta", "Pseudomonadota", "Bacteroidota"]
    label_by_contig = {
        "ctg1_host_A": "Chlorophyta",
        "ctg2_host_A_hap": "Chlorophyta",
        "ctg3_host_B": "Chlorophyta",
        "ctg6_host_C": "Chlorophyta",
        "ctg4_contam_bac1": "Pseudomonadota",
        "ctg5_contam_bac2": "Bacteroidota",
        # Deliberately the same excluded taxon as ctg4_contam_bac1 -- see the
        # module docstring above. This is what makes the fixture a regression
        # test for ORGANELLE_ISOLATE running before BTK_FILTER: mitochondria
        # are frequently taxonomically misclassified as bacterial
        # contamination in real blob plots, and this contig must survive
        # BTK_FILTER's taxonomic filter regardless of this label.
        "ctg7_organelle": "Pseudomonadota",
    }
    values = [keys.index(label_by_contig[i]) for i in identifiers]

    # Compression is deliberately mixed here (see write_blobdir_json): a real
    # sanger-tol/blobtoolkit BlobDir gzip-compresses its per-field JSON
    # files, but not every BlobDir out there will be -- mixing both in this
    # one fixture exercises both code paths without a second nf-test scenario.
    write_blobdir_json(blobdir / "meta.json", {
        "id": "blobpurge_test",
        "name": "blobpurge synthetic test BlobDir",
        "record_type": "contig",
        "settings": {"taxrule": "bestsumorder", "taxrules": ["bestsumorder"]},
        "fields": [
            {"id": "identifiers", "type": "identifier"},
            {"id": "length", "type": "variable"},
            {"id": "gc", "type": "variable"},
            {"id": "ncount", "type": "variable"},
            {"id": "buscoregions_phylum", "type": "category"},
            # blobtools filter's own --summary STDOUT overview always looks up
            # a "<taxrule>_<summary-rank>" category field (default rank:
            # phylum), independently of whatever field --param filters on.
            # Real sanger-tol/blobtoolkit BlobDirs always have this from their
            # own taxonomy-hit summarisation step.
            {"id": "bestsumorder_phylum", "type": "category"},
        ],
    }, gzip_compress=True)

    write_blobdir_json(blobdir / "identifiers.json", {
        "id": "identifiers", "name": "identifiers", "type": "identifier", "values": identifiers,
    }, gzip_compress=False)

    write_blobdir_json(blobdir / "length.json", {
        "id": "length", "name": "length", "type": "variable", "values": lengths,
    }, gzip_compress=True)

    write_blobdir_json(blobdir / "gc.json", {
        "id": "gc", "name": "GC", "type": "variable", "values": gc_portions,
    }, gzip_compress=False)

    write_blobdir_json(blobdir / "ncount.json", {
        "id": "ncount", "name": "N count", "type": "variable", "values": n_counts,
    }, gzip_compress=False)

    write_blobdir_json(blobdir / "buscoregions_phylum.json", {
        "id": "buscoregions_phylum", "name": "buscoregions_phylum", "type": "category",
        "keys": keys, "values": values,
    }, gzip_compress=True)

    write_blobdir_json(blobdir / "bestsumorder_phylum.json", {
        "id": "bestsumorder_phylum", "name": "bestsumorder_phylum", "type": "category",
        "keys": keys, "values": values,
    }, gzip_compress=False)

    # No samplesheet.csv is written here: nf-schema resolves relative paths
    # in a samplesheet against the launch directory, not the CSV's own
    # location, so a portable samplesheet needs absolute paths. conf/test.config
    # builds one at config-load time from ${projectDir}; see there.

    print(f"Wrote {len(records)} contigs, {len(all_pairs)} read pairs.")


if __name__ == "__main__":
    main()
