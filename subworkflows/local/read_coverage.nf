//
// READ_COVERAGE: produce a purge_dups (ngscstat) coverage stat/base-cov pair for each sample.
//
// reads_cram is already aligned upstream to the *original*, unfiltered
// assembly (the sanger-tol/blobtoolkit precondition for building the
// BlobDir's coverage tracks); it is subset down to the contigs retained by
// BTK_FILTER. ngscstat is then run on the resulting sorted/indexed BAM --
// purge_dups' Illumina coverage path, which is explicitly logged as less
// validated than the PacBio path (see PURGEDUPS_NGSCSTAT).
//

include { CRAM_TO_BAM                         } from '../../modules/local/cram_to_bam/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME } from '../../modules/nf-core/samtools/sort/main'
include { BAM_STATS_SAMTOOLS                  } from '../../subworkflows/nf-core/bam_stats_samtools/main'
include { PURGEDUPS_NGSCSTAT                  } from '../../modules/local/purgedups_ngscstat/main'

workflow READ_COVERAGE {

    take:
    ch_reads
    // tuple( meta, original_assembly, reads_cram, retained_ids )

    main:
    ch_versions = Channel.empty()
    ch_caveats  = Channel.empty()

    ch_no_fasta_fai = Channel.value([ [:], [], [] ])

    CRAM_TO_BAM(ch_reads)
    ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)

    ch_bam = CRAM_TO_BAM.out.bam

    BAM_STATS_SAMTOOLS(ch_bam, ch_no_fasta_fai)

    //
    // ngscstat needs a name-sorted BAM (see SAMTOOLS_SORT_NAME); the
    // coordinate-sorted/indexed BAM is kept as-is for downstream consumers
    // (e.g. purge_haplotigs).
    //
    SAMTOOLS_SORT_NAME(ch_bam.map { meta, bam, _bai -> [ meta, bam ] }, ch_no_fasta_fai, '')

    PURGEDUPS_NGSCSTAT(SAMTOOLS_SORT_NAME.out.bam)
    ch_versions = ch_versions.mix(PURGEDUPS_NGSCSTAT.out.versions)
    ch_caveats  = ch_caveats.mix(PURGEDUPS_NGSCSTAT.out.caveat)

    //
    // QC-only: mapping-rate/insert-size/per-contig stats for MultiQC.
    //
    ch_multiqc_files = BAM_STATS_SAMTOOLS.out.stats.map { _meta, f -> f }
        .mix(BAM_STATS_SAMTOOLS.out.flagstat.map { _meta, f -> f })
        .mix(BAM_STATS_SAMTOOLS.out.idxstats.map { _meta, f -> f })

    emit:
    stat          = PURGEDUPS_NGSCSTAT.out.stat
    base_cov      = PURGEDUPS_NGSCSTAT.out.base_cov
    bam           = ch_bam
    caveats       = ch_caveats
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
