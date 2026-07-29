#!/usr/bin/env python3
"""Generate the tiny synthetic dataset used by -profile test.

Produces:
  - assembly.fasta.gz: 9 contigs -- two heterozygous duplicate pairs
    (ctg1/ctg2 and ctg7/ctg8, each 8% divergent, so BTK_FILTER->purge_dups
    and HaploMerger2 both have something real to purge/self-align), three
    unrelated unique host contigs, and two "contaminant" contigs. Two
    duplicate pairs (rather than one) give HaploMerger2's self-alignment
    graph enough scaffolds/aligning pairs post-BTK_FILTER for its own
    HM_pathFinder_preparation.pl to build a workable graph -- a single pair
    (4 scaffolds, 1 aligning pair) was too sparse and crashed it.
  - blobdir/: a minimal, hand-built BlobDir (meta.json + identifiers/length/
    buscoregions_phylum field JSON) matching the schema BTK_FILTER expects.
  - reads_R1.fastq.gz / reads_R2.fastq.gz: paired-end reads sampled
    independently from each half of both duplicate pairs (ctg1/ctg2 and
    ctg7/ctg8, each at half the depth used for the unique host contigs) --
    reflecting the real biology of an uncollapsed
    heterozygous locus, where each assembled copy only receives reads from
    the haplotype it truly represents. At 8% divergence each read maps
    confidently (MAPQ > 30) to its true contig of origin -- required for
    purge_dups' ngscstat, which discards low-MAPQ reads -- while remaining
    well within minimap2 asm5's detection range for the self-alignment step.
    An earlier design tried near-zero divergence to force bwa-mem2
    multi-mapping ambiguity, but that MAPQ~0 signal is exactly what
    ngscstat's default -q 30 filter discards, silently zeroing out the
    coverage it was meant to produce. No reads are simulated for the
    contaminant contigs.
  - samplesheet.csv

Everything here is synthetic and deterministic (fixed RNG seed) -- this is a
plumbing fixture for -profile test, not a biologically meaningful assembly.
"""
import gzip
import json
import random
from pathlib import Path

random.seed(42)
OUT = Path(__file__).parent
BASES = "ACGT"


def random_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


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
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for rid, seq in records.items():
            handle.write(f">{rid}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")


def sample_reads(seq, contig_id, depth, read_len=100, insert=300, error_rate=0.005):
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
        r1 = frag[:read_len]
        r2 = revcomp(frag[-read_len:])
        r1 = mutate(r1, error_rate)
        r2 = mutate(r2, error_rate)
        pairs.append((f"{contig_id}_{i}", r1, r2))
    return pairs


def write_fastq(path, reads, read_index):
    with gzip.open(path, "wt") as handle:
        for name, seq in reads:
            qual = "I" * len(seq)
            handle.write(f"@{name}/{read_index}\n{seq}\n+\n{qual}\n")


def main():
    ctg1 = random_seq(3000)
    ctg2 = mutate(ctg1, 0.08)
    ctg7 = random_seq(3000)
    ctg8 = mutate(ctg7, 0.08)
    # Unique host contigs are kept substantially larger than the duplicated
    # loci so the assembly-wide depth histogram's tallest peak is the
    # normal single-copy coverage level, not the half-depth haplotig region
    # -- matching what purge_haplotigs' cutoff estimation assumes (most of a
    # real assembly is single-copy; only a minority is uncollapsed haplotig).
    # Two duplicate pairs (rather than one) also give HaploMerger2's
    # self-alignment graph enough scaffolds/aligning pairs post-BTK_FILTER
    # for HM_pathFinder_preparation.pl to build a workable graph.
    ctg3 = random_seq(12000)
    ctg6 = random_seq(8000)
    ctg9 = random_seq(20000)
    ctg4 = random_seq(1200)
    ctg5 = random_seq(1000)

    records = {
        "ctg1_host_A": ctg1,
        "ctg2_host_A_hap": ctg2,
        "ctg3_host_B": ctg3,
        "ctg4_contam_bac1": ctg4,
        "ctg5_contam_bac2": ctg5,
        "ctg6_host_C": ctg6,
        "ctg7_host_D": ctg7,
        "ctg8_host_D_hap": ctg8,
        "ctg9_host_E": ctg9,
    }
    write_fasta(OUT / "assembly.fasta.gz", records)

    depth = 80
    all_pairs = []
    all_pairs += sample_reads(ctg1, "h1", depth=depth // 2)
    all_pairs += sample_reads(ctg2, "h2", depth=depth // 2)
    all_pairs += sample_reads(ctg3, "ctg3", depth=depth)
    all_pairs += sample_reads(ctg6, "ctg6", depth=depth)
    all_pairs += sample_reads(ctg7, "h3", depth=depth // 2)
    all_pairs += sample_reads(ctg8, "h4", depth=depth // 2)
    all_pairs += sample_reads(ctg9, "ctg9", depth=depth)

    r1_reads = [(name, r1) for name, r1, _r2 in all_pairs]
    r2_reads = [(name, r2) for name, _r1, r2 in all_pairs]
    write_fastq(OUT / "reads_R1.fastq.gz", r1_reads, 1)
    write_fastq(OUT / "reads_R2.fastq.gz", r2_reads, 2)

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
        "ctg7_host_D": "Chlorophyta",
        "ctg8_host_D_hap": "Chlorophyta",
        "ctg9_host_E": "Chlorophyta",
        "ctg4_contam_bac1": "Pseudomonadota",
        "ctg5_contam_bac2": "Bacteroidota",
    }
    values = [keys.index(label_by_contig[i]) for i in identifiers]

    (blobdir / "meta.json").write_text(json.dumps({
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
    }, indent=2))

    (blobdir / "identifiers.json").write_text(json.dumps({
        "id": "identifiers", "name": "identifiers", "type": "identifier", "values": identifiers,
    }, indent=2))

    (blobdir / "length.json").write_text(json.dumps({
        "id": "length", "name": "length", "type": "variable", "values": lengths,
    }, indent=2))

    (blobdir / "gc.json").write_text(json.dumps({
        "id": "gc", "name": "GC", "type": "variable", "values": gc_portions,
    }, indent=2))

    (blobdir / "ncount.json").write_text(json.dumps({
        "id": "ncount", "name": "N count", "type": "variable", "values": n_counts,
    }, indent=2))

    (blobdir / "buscoregions_phylum.json").write_text(json.dumps({
        "id": "buscoregions_phylum", "name": "buscoregions_phylum", "type": "category",
        "keys": keys, "values": values,
    }, indent=2))

    (blobdir / "bestsumorder_phylum.json").write_text(json.dumps({
        "id": "bestsumorder_phylum", "name": "bestsumorder_phylum", "type": "category",
        "keys": keys, "values": values,
    }, indent=2))

    # No samplesheet.csv is written here: nf-schema resolves relative paths
    # in a samplesheet against the launch directory, not the CSV's own
    # location, so a portable samplesheet needs absolute paths. conf/test.config
    # builds one at config-load time from ${projectDir}; see there.

    print(f"Wrote {len(records)} contigs, {len(all_pairs)} read pairs.")


if __name__ == "__main__":
    main()
