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

include { CRAM_TO_BAM         } from '../../modules/local/cram_to_bam/main'
include { BWAMEM2_INDEX       } from '../../modules/local/bwamem2_index/main'
include { BWAMEM2_MEM         } from '../../modules/local/bwamem2_mem/main'
include { SAMTOOLS_SORT_INDEX } from '../../modules/local/samtools_sort_index/main'
include { SAMTOOLS_SORT_NAME  } from '../../modules/local/samtools_sort_name/main'
include { SAMTOOLS_STATS      } from '../../modules/local/samtools_stats/main'
include { PURGEDUPS_NGSCSTAT  } from '../../modules/local/purgedups_ngscstat/main'

workflow READ_COVERAGE {

    take:
    ch_reads
    // tuple( meta, original_assembly, filtered_assembly, retained_ids, reads_cram, reads_r1, reads_r2 )

    main:
    ch_versions = Channel.empty()
    ch_caveats  = Channel.empty()

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

    //
    // FASTQ path: fresh bwa-mem2 alignment against the filtered assembly.
    //
    ch_branched.fastq
        .map { meta, filtered_assembly, reads_r1, reads_r2 -> [ meta, filtered_assembly ] }
        .set { ch_for_index }

    BWAMEM2_INDEX(ch_for_index)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    ch_branched.fastq
        .join(BWAMEM2_INDEX.out.index)
        .map { meta, filtered_assembly, reads_r1, reads_r2, index -> [ meta, index, reads_r1, reads_r2 ] }
        .set { ch_for_mem }

    BWAMEM2_MEM(ch_for_mem)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    SAMTOOLS_SORT_INDEX(BWAMEM2_MEM.out.sam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions)

    //
    // Common: ngscstat on whichever BAM was produced. ngscstat needs a
    // name-sorted BAM (see SAMTOOLS_SORT_NAME); the coordinate-sorted/indexed
    // BAM is kept as-is for downstream consumers (e.g. purge_haplotigs).
    //
    ch_bam = CRAM_TO_BAM.out.bam.mix(SAMTOOLS_SORT_INDEX.out.bam)

    SAMTOOLS_SORT_NAME(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_NAME.out.versions)

    PURGEDUPS_NGSCSTAT(SAMTOOLS_SORT_NAME.out.bam)
    ch_versions = ch_versions.mix(PURGEDUPS_NGSCSTAT.out.versions)
    ch_caveats  = ch_caveats.mix(PURGEDUPS_NGSCSTAT.out.caveat)

    //
    // QC-only: mapping-rate/insert-size/per-contig stats for MultiQC, on
    // whichever BAM was produced (CRAM subset or fresh bwa-mem2 mapping).
    //
    SAMTOOLS_STATS(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    ch_multiqc_files = SAMTOOLS_STATS.out.stats.map { _meta, f -> f }
        .mix(SAMTOOLS_STATS.out.flagstat.map { _meta, f -> f })
        .mix(SAMTOOLS_STATS.out.idxstats.map { _meta, f -> f })

    emit:
    stat          = PURGEDUPS_NGSCSTAT.out.stat
    base_cov      = PURGEDUPS_NGSCSTAT.out.base_cov
    bam           = ch_bam
    caveats       = ch_caveats
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
