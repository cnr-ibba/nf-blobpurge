//
// READ_COVERAGE: produce a purge_dups (ngscstat) coverage stat/base-cov pair for each sample.
//
// If reads_cram is available, the CRAM (already aligned upstream to the
// *original*, unfiltered assembly) is subset down to the contigs retained by
// BTK_FILTER. Otherwise, reads_r1/reads_r2 are mapped fresh onto the
// *filtered* assembly with bwa-mem2. Either way, ngscstat is run on the
// resulting sorted/indexed BAM -- purge_dups' Illumina coverage path, which
// is explicitly logged as less validated than the PacBio path (see
// PURGEDUPS_NGSCSTAT).
//

include { CRAM_TO_BAM                         } from '../../modules/local/cram_to_bam/main'
include { BWAMEM2_INDEX                       } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM                         } from '../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME } from '../../modules/nf-core/samtools/sort/main'
include { BAM_SORT_STATS_SAMTOOLS             } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_STATS_SAMTOOLS                  } from '../../subworkflows/nf-core/bam_stats_samtools/main'
include { PURGEDUPS_NGSCSTAT                  } from '../../modules/local/purgedups_ngscstat/main'

workflow READ_COVERAGE {

    take:
    ch_reads
    // tuple( meta, original_assembly, filtered_assembly, retained_ids, reads_cram, reads_r1, reads_r2 )

    main:
    ch_versions = Channel.empty()
    ch_caveats  = Channel.empty()

    // No fasta reference needed by any of the nf-core samtools/bwamem2 calls
    // below (all consume/produce BAM, never CRAM), so a constant empty value
    // channel is broadcast to every one of them.
    ch_no_fasta     = Channel.value([ [:], [] ])
    ch_no_fasta_fai = Channel.value([ [:], [], [] ])

    ch_reads
        .branch { meta, original_assembly, filtered_assembly, retained_ids, reads_cram, reads_r1, reads_r2 ->
            cram:  meta.has_cram
                return [ meta, original_assembly, reads_cram, retained_ids ]
            fastq: !meta.has_cram
                return [ meta, filtered_assembly, reads_r1, reads_r2 ]
        }
        .set { ch_branched }

    //
    // CRAM path: subset the existing alignment, no remapping.
    //
    CRAM_TO_BAM(ch_branched.cram)
    ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)

    BAM_STATS_SAMTOOLS(CRAM_TO_BAM.out.bam, ch_no_fasta_fai)

    //
    // FASTQ path: fresh bwa-mem2 alignment against the filtered assembly.
    //
    ch_branched.fastq
        .map { meta, filtered_assembly, reads_r1, reads_r2 -> [ meta, filtered_assembly ] }
        .set { ch_for_index }

    BWAMEM2_INDEX(ch_for_index)

    ch_branched.fastq
        .join(BWAMEM2_INDEX.out.index)
        .multiMap { meta, filtered_assembly, reads_r1, reads_r2, index ->
            reads: [ meta, [ reads_r1, reads_r2 ] ]
            index: [ meta, index ]
        }
        .set { ch_for_mem }

    BWAMEM2_MEM(ch_for_mem.reads, ch_for_mem.index, ch_no_fasta, false)

    BAM_SORT_STATS_SAMTOOLS(BWAMEM2_MEM.out.bam, ch_no_fasta_fai)

    //
    // Common: ngscstat on whichever BAM was produced. ngscstat needs a
    // name-sorted BAM (see SAMTOOLS_SORT_NAME); the coordinate-sorted/indexed
    // BAM is kept as-is for downstream consumers (e.g. purge_haplotigs).
    //
    ch_bam = CRAM_TO_BAM.out.bam.mix(
        BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.index)
    )

    SAMTOOLS_SORT_NAME(ch_bam.map { meta, bam, _bai -> [ meta, bam ] }, ch_no_fasta_fai, '')

    PURGEDUPS_NGSCSTAT(SAMTOOLS_SORT_NAME.out.bam)
    ch_versions = ch_versions.mix(PURGEDUPS_NGSCSTAT.out.versions)
    ch_caveats  = ch_caveats.mix(PURGEDUPS_NGSCSTAT.out.caveat)

    //
    // QC-only: mapping-rate/insert-size/per-contig stats for MultiQC, on
    // whichever BAM was produced (CRAM subset or fresh bwa-mem2 mapping).
    //
    ch_multiqc_files = BAM_STATS_SAMTOOLS.out.stats.map { _meta, f -> f }
        .mix(BAM_STATS_SAMTOOLS.out.flagstat.map { _meta, f -> f })
        .mix(BAM_STATS_SAMTOOLS.out.idxstats.map { _meta, f -> f })
        .mix(BAM_SORT_STATS_SAMTOOLS.out.stats.map { _meta, f -> f })
        .mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.map { _meta, f -> f })
        .mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.map { _meta, f -> f })

    emit:
    stat          = PURGEDUPS_NGSCSTAT.out.stat
    base_cov      = PURGEDUPS_NGSCSTAT.out.base_cov
    bam           = ch_bam
    caveats       = ch_caveats
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
