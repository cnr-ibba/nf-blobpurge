//
// ORGANELLE_ISOLATE: isolate organelle-like (mitochondrial/plastid) contigs,
// by coverage (+ optional reference-match corroboration), from the RAW
// assembly, before BTK_FILTER ever runs.
//
// Mitochondria are frequently taxonomically misclassified as bacterial
// contamination in blob plots, so running this classification *after*
// BTK_FILTER (an earlier design of this subworkflow) cannot protect an
// organelle contig that blobtools filter already discarded as a
// "contaminant" -- by the time this subworkflow would see it, it's gone.
// Classifying on the raw assembly instead, using coverage from the reads
// CRAM that's already aligned to it (sanger-tol/blobtoolkit's own
// precondition -- no re-mapping needed), and physically splitting organelle
// contigs out before BTK_FILTER's `--fasta` input is built, means
// blobtools' taxonomic filter can never see -- and so can never discard --
// an organelle contig.
//
// Neither purge_dups' .stat file nor purge_haplotigs' samtools-depth TSV has
// any per-contig awareness -- both are genome-wide pooled histograms -- so
// keeping organelle contigs (10-100x nuclear coverage) out of BTK_FILTER's
// (and therefore READ_COVERAGE's) output is also what keeps them out of
// calcuts/estimate_purgehaplotigs_cutoffs.py's cutoff estimation, with no
// separate accounting needed downstream.
//
// Isolated contigs are not contamination: they are re-merged into each
// stage's FASTA downstream (see workflows/nf-blobpurge.nf), and a matching
// paired-end read subset is exported for use with external organelle-assembly
// tools (GetOrganelle, MitoHiFi, oatk, ...) -- this pipeline does not attempt
// to assemble organelles itself.
//

include { CRAM_TO_BAM as CRAM_TO_BAM_RAW                   } from '../../modules/local/cram_to_bam/main'
include { SAMTOOLS_COVERAGE                                } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME_ORGANELLE     } from '../../modules/nf-core/samtools/sort/main'
include { ORGANELLE_REFERENCE_ALIGN                         } from '../../modules/local/organelle_reference_align/main'
include { ORGANELLE_DETECT                                  } from '../../modules/local/organelle_detect/main'
include { SAMTOOLS_VIEW_SUBSET as ORGANELLE_FILTER_BAM      } from '../../modules/local/samtools_view_subset/main'
include { ORGANELLE_EXTRACT_READS                           } from '../../modules/local/organelle_extract_reads/main'

workflow ORGANELLE_ISOLATE {

    take:
    ch_reads     // tuple( meta, raw_assembly, reads_cram )
    ch_reference // tuple( meta, reference_fasta_or_NO_FILE )

    main:
    ch_versions       = Channel.empty()
    ch_caveats        = Channel.empty()
    ch_multiqc_files  = Channel.empty()
    ch_no_fasta_fai   = Channel.value([ [:], [], [] ])

    //
    // Raw-assembly-wide coverage BAM. CRAM_TO_BAM_RAW is the same module
    // READ_COVERAGE uses to subset the CRAM onto BTK_FILTER-retained
    // contigs, aliased here and fed assets/NO_FILE (a 0-byte sentinel) as
    // its "retained_ids" argument: with an empty ID list its `samtools view`
    // call gets zero region arguments, which returns every contig -- i.e.
    // "no filtering" -- so no separate module is needed to get a raw,
    // whole-assembly BAM out of the same CRAM.
    //
    ch_reads
        .combine(Channel.fromPath("${projectDir}/assets/NO_FILE"))
        .set { ch_for_raw_bam }

    CRAM_TO_BAM_RAW(ch_for_raw_bam)
    ch_versions = ch_versions.mix(CRAM_TO_BAM_RAW.out.versions)

    //
    // Per-contig mean depth on the raw assembly -- no filtering decision has
    // been made yet at this point, so every contig is scored.
    //
    SAMTOOLS_COVERAGE(CRAM_TO_BAM_RAW.out.bam, ch_no_fasta_fai)

    //
    // Optional corroborating evidence: align the raw (still-unclassified)
    // assembly against a user-supplied organelle reference. Skipped (empty
    // PAF) when no reference is given.
    //
    def ch_raw_fasta = ch_reads.map { meta, assembly, _cram -> [ meta, assembly ] }

    ORGANELLE_REFERENCE_ALIGN(ch_raw_fasta.join(ch_reference))
    ch_versions = ch_versions.mix(ORGANELLE_REFERENCE_ALIGN.out.versions)

    //
    // Classification + FASTA split.
    //
    ORGANELLE_DETECT(
        ch_raw_fasta
            .join(SAMTOOLS_COVERAGE.out.coverage)
            .join(ORGANELLE_REFERENCE_ALIGN.out.paf)
    )
    ch_versions = ch_versions.mix(ORGANELLE_DETECT.out.versions)
    ch_caveats  = ch_caveats.mix(ORGANELLE_DETECT.out.caveat)
    ch_multiqc_files = ch_multiqc_files.mix(ORGANELLE_DETECT.out.mqc_json.map { _meta, f -> f })

    //
    // Export the paired-end reads mapping to the isolated organelle contigs.
    // (The nuclear side needs no equivalent BAM split here: BTK_FILTER now
    // runs on ORGANELLE_DETECT.out.nuclear_fasta directly, so READ_COVERAGE
    // downstream is automatically organelle-free by construction.)
    //
    ORGANELLE_FILTER_BAM(CRAM_TO_BAM_RAW.out.bam.join(ORGANELLE_DETECT.out.organelle_ids))
    ch_versions = ch_versions.mix(ORGANELLE_FILTER_BAM.out.versions)

    SAMTOOLS_SORT_NAME_ORGANELLE(
        ORGANELLE_FILTER_BAM.out.bam.map { meta, bam, _bai -> [ meta, bam ] },
        ch_no_fasta_fai,
        ''
    )
    ORGANELLE_EXTRACT_READS(SAMTOOLS_SORT_NAME_ORGANELLE.out.bam)
    ch_versions = ch_versions.mix(ORGANELLE_EXTRACT_READS.out.versions)

    emit:
    nuclear_fasta   = ORGANELLE_DETECT.out.nuclear_fasta    // tuple(meta, fasta) -> BTK_FILTER's `assembly` input
    organelle_fasta = ORGANELLE_DETECT.out.organelle_fasta  // tuple(meta, fasta) -> re-merged into each stage's FASTA
    organelle_ids   = ORGANELLE_DETECT.out.organelle_ids    // tuple(meta, ids)   -> BTK_FILTER's check_span --exclude-ids
    organelle_reads = ORGANELLE_EXTRACT_READS.out.reads     // tuple(meta, r1, r2) -> published output only, not consumed downstream
    report_json     = ORGANELLE_DETECT.out.report_json      // tuple(meta, json) -> BLOBPURGE_REPORT input
    caveats         = ch_caveats
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
