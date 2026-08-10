//
// ORGANELLE_ISOLATE: isolate organelle-like (mitochondrial/plastid) contigs,
// by coverage (+ optional reference-match corroboration), before they reach
// PURGE_DUPS/PURGE_HAPLOTIGS.
//
// Neither purge_dups' .stat file nor purge_haplotigs' samtools-depth TSV has
// any per-contig awareness -- both are genome-wide pooled histograms -- so
// the only clean way to keep organelle contigs (10-100x nuclear coverage)
// from distorting calcuts/estimate_purgehaplotigs_cutoffs.py's cutoffs is to
// keep their positions out of those files in the first place, not to try to
// "clean" the pooled histograms after the fact. READ_COVERAGE itself is left
// untouched (re-mapping/re-subsetting the CRAM just to exclude a handful of
// contigs would be wasted, expensive work): this subworkflow instead works
// from READ_COVERAGE's own BAM, classifies contigs, filters the BAM down to
// nuclear-only, and re-runs PURGEDUPS_NGSCSTAT (unchanged) on that filtered,
// name-sorted BAM to get a genuinely organelle-free .stat/base_cov pair.
//
// Isolated contigs are not contamination: they are re-merged into each purge
// stage's FASTA downstream (see workflows/nf-blobpurge.nf), and a matching
// paired-end read subset is exported for use with external organelle-assembly
// tools (GetOrganelle, MitoHiFi, oatk, ...) -- this pipeline does not attempt
// to assemble organelles itself.
//

include { SAMTOOLS_COVERAGE                                } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME_NUCLEAR       } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME_ORGANELLE     } from '../../modules/nf-core/samtools/sort/main'
include { PURGEDUPS_NGSCSTAT                                } from '../../modules/local/purgedups_ngscstat/main'
include { ORGANELLE_REFERENCE_ALIGN                         } from '../../modules/local/organelle_reference_align/main'
include { ORGANELLE_DETECT                                  } from '../../modules/local/organelle_detect/main'
include { SAMTOOLS_VIEW_SUBSET as NUCLEAR_FILTER_BAM        } from '../../modules/local/samtools_view_subset/main'
include { SAMTOOLS_VIEW_SUBSET as ORGANELLE_FILTER_BAM      } from '../../modules/local/samtools_view_subset/main'
include { ORGANELLE_EXTRACT_READS                           } from '../../modules/local/organelle_extract_reads/main'

workflow ORGANELLE_ISOLATE {

    take:
    ch_fasta     // tuple( meta, filtered_fasta )              -- BTK_FILTER.out.fasta
    ch_bam       // tuple( meta, bam, bai )                     -- READ_COVERAGE.out.bam (coord-sorted+indexed)
    ch_reference // tuple( meta, reference_fasta_or_NO_FILE )

    main:
    ch_versions       = Channel.empty()
    ch_caveats        = Channel.empty()
    ch_multiqc_files  = Channel.empty()
    ch_no_fasta_fai   = Channel.value([ [:], [], [] ])

    //
    // Per-contig mean depth from the BAM READ_COVERAGE already produced --
    // no re-mapping needed to make a classification decision.
    //
    SAMTOOLS_COVERAGE(ch_bam, ch_no_fasta_fai)

    //
    // Optional corroborating evidence: align the (still-unclassified)
    // filtered assembly against a user-supplied organelle reference.
    // Skipped (empty PAF) when no reference is given.
    //
    ORGANELLE_REFERENCE_ALIGN(ch_fasta.join(ch_reference))
    ch_versions = ch_versions.mix(ORGANELLE_REFERENCE_ALIGN.out.versions)

    //
    // Classification + FASTA split.
    //
    ORGANELLE_DETECT(
        ch_fasta
            .join(SAMTOOLS_COVERAGE.out.coverage)
            .join(ORGANELLE_REFERENCE_ALIGN.out.paf)
    )
    ch_versions = ch_versions.mix(ORGANELLE_DETECT.out.versions)
    ch_caveats  = ch_caveats.mix(ORGANELLE_DETECT.out.caveat)
    ch_multiqc_files = ch_multiqc_files.mix(ORGANELLE_DETECT.out.mqc_json.map { _meta, f -> f })

    //
    // Subset the coordinate-sorted BAM by contig-ID list, both ways: nuclear
    // contigs feed ngscstat (below) and PURGE_HAPLOTIGS; organelle contigs
    // feed the read-export step.
    //
    NUCLEAR_FILTER_BAM(ch_bam.join(ORGANELLE_DETECT.out.nuclear_ids))
    ch_versions = ch_versions.mix(NUCLEAR_FILTER_BAM.out.versions)

    ORGANELLE_FILTER_BAM(ch_bam.join(ORGANELLE_DETECT.out.organelle_ids))
    ch_versions = ch_versions.mix(ORGANELLE_FILTER_BAM.out.versions)

    //
    // Re-derive an organelle-free .stat/base_cov pair: ngscstat (like
    // calcuts/hist_plot.py downstream) needs a name-sorted BAM, not the
    // coordinate-sorted one used everywhere else -- see PURGEDUPS_NGSCSTAT.
    // The Illumina-coverage-path caveat was already surfaced once by
    // READ_COVERAGE's own ngscstat run, so it is deliberately not re-emitted
    // here to avoid an identical duplicate bullet in the report.
    //
    SAMTOOLS_SORT_NAME_NUCLEAR(
        NUCLEAR_FILTER_BAM.out.bam.map { meta, bam, _bai -> [ meta, bam ] },
        ch_no_fasta_fai,
        ''
    )
    PURGEDUPS_NGSCSTAT(SAMTOOLS_SORT_NAME_NUCLEAR.out.bam)
    ch_versions = ch_versions.mix(PURGEDUPS_NGSCSTAT.out.versions)

    //
    // Export the paired-end reads mapping to the isolated organelle contigs.
    //
    SAMTOOLS_SORT_NAME_ORGANELLE(
        ORGANELLE_FILTER_BAM.out.bam.map { meta, bam, _bai -> [ meta, bam ] },
        ch_no_fasta_fai,
        ''
    )
    ORGANELLE_EXTRACT_READS(SAMTOOLS_SORT_NAME_ORGANELLE.out.bam)
    ch_versions = ch_versions.mix(ORGANELLE_EXTRACT_READS.out.versions)

    emit:
    nuclear_fasta   = ORGANELLE_DETECT.out.nuclear_fasta    // tuple(meta, fasta) -> replaces BTK_FILTER.out.fasta into PURGE_DUPS/PURGE_HAPLOTIGS
    organelle_fasta = ORGANELLE_DETECT.out.organelle_fasta  // tuple(meta, fasta) -> re-merged into each purge stage's FASTA
    stat            = PURGEDUPS_NGSCSTAT.out.stat           // tuple(meta, stat) -> replaces READ_COVERAGE.out.stat into PURGE_DUPS
    base_cov        = PURGEDUPS_NGSCSTAT.out.base_cov       // tuple(meta, base_cov) -> replaces READ_COVERAGE.out.base_cov into PURGE_DUPS
    bam             = NUCLEAR_FILTER_BAM.out.bam             // tuple(meta, bam, bai) -> replaces READ_COVERAGE.out.bam into PURGE_HAPLOTIGS
    organelle_reads = ORGANELLE_EXTRACT_READS.out.reads      // tuple(meta, r1, r2) -> published output only, not consumed downstream
    report_json     = ORGANELLE_DETECT.out.report_json       // tuple(meta, json) -> new BLOBPURGE_REPORT input
    caveats         = ch_caveats
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
