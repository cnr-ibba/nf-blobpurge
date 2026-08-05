//
// PURGE_HAPLOTIGS: independent cross-check of PURGE_DUPS.
//
// Cutoffs for `purge_haplotigs cov` are normally chosen by eye from the
// coverage histogram; here they are estimated automatically from the
// sample's own depth distribution (see PURGEHAPLOTIGS_ESTIMATE_CUTOFFS).
// The coverage BAM reused here comes from READ_COVERAGE, which for most
// samples will be short-read (Illumina) coverage -- purge_haplotigs is
// designed and documented primarily around long-read coverage, so this is
// logged as a caveat rather than presented as an equivalent use case.
//

include { SAMTOOLS_DEPTH                   } from '../../modules/local/samtools_depth/main'
include { PURGEHAPLOTIGS_ESTIMATE_CUTOFFS  } from '../../modules/local/purgehaplotigs_estimate_cutoffs/main'
include { PURGEHAPLOTIGS_HIST              } from '../../modules/local/purgehaplotigs_hist/main'
include { PURGEHAPLOTIGS_COV                } from '../../modules/local/purgehaplotigs_cov/main'
include { PURGEHAPLOTIGS_PURGE              } from '../../modules/local/purgehaplotigs_purge/main'

workflow PURGE_HAPLOTIGS {

    take:
    ch_fasta // tuple( meta, filtered_fasta )
    ch_bam   // tuple( meta, bam, bai )

    main:
    ch_versions = Channel.empty()
    ch_caveats  = Channel.empty()

    ch_caveats = ch_caveats.mix(
        ch_bam.map { meta, bam, bai ->
            [ meta, "purge_haplotigs cross-check for '${meta.id}' was run on the same coverage used for purge_dups (typically short-read Illumina): purge_haplotigs is designed and documented primarily for long-read coverage, so this should be considered an indicative cross-check, not equivalent to its standard use case." ]
        }
    )

    SAMTOOLS_DEPTH(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    PURGEHAPLOTIGS_ESTIMATE_CUTOFFS(SAMTOOLS_DEPTH.out.depth)
    ch_versions = ch_versions.mix(PURGEHAPLOTIGS_ESTIMATE_CUTOFFS.out.versions)

    ch_bam
        .join(ch_fasta)
        .set { ch_for_hist }

    PURGEHAPLOTIGS_HIST(ch_for_hist)
    ch_versions = ch_versions.mix(PURGEHAPLOTIGS_HIST.out.versions)

    PURGEHAPLOTIGS_HIST.out.gencov
        .join(PURGEHAPLOTIGS_ESTIMATE_CUTOFFS.out.cutoffs)
        .map { meta, gencov, cutoffs_json ->
            def cutoffs = new groovy.json.JsonSlurper().parse(cutoffs_json)
            [ meta, gencov, cutoffs.low, cutoffs.mid, cutoffs.high ]
        }
        .set { ch_for_cov }

    PURGEHAPLOTIGS_COV(ch_for_cov)
    ch_versions = ch_versions.mix(PURGEHAPLOTIGS_COV.out.versions)

    ch_fasta
        .join(PURGEHAPLOTIGS_COV.out.coverage_stats)
        .set { ch_for_purge }

    PURGEHAPLOTIGS_PURGE(ch_for_purge)
    ch_versions = ch_versions.mix(PURGEHAPLOTIGS_PURGE.out.versions)

    emit:
    purged_fasta  = PURGEHAPLOTIGS_PURGE.out.purged_fasta
    removed_fasta = PURGEHAPLOTIGS_PURGE.out.removed_fasta
    histogram_png = PURGEHAPLOTIGS_HIST.out.histogram_png
    cutoffs       = PURGEHAPLOTIGS_ESTIMATE_CUTOFFS.out.cutoffs
    caveats       = ch_caveats
    versions      = ch_versions
}
