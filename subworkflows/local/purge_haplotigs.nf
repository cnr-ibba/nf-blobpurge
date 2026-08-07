//
// PURGE_HAPLOTIGS: independent cross-check of PURGE_DUPS.
//
// Cutoffs for `purge_haplotigs cov` are normally chosen by eye from the
// coverage histogram; here they are estimated automatically from the
// sample's own depth distribution (see PURGEHAPLOTIGS_ESTIMATE_CUTOFFS). If
// that estimation finds no usable bimodal signal for a sample, it reports
// itself as skipped rather than guessing, and this subworkflow simply omits
// that sample's 'purge_haplotigs' stage (with an explanatory caveat) instead
// of failing the whole run.
//
// The coverage BAM reused here comes from READ_COVERAGE, which for most
// samples will be short-read (Illumina) coverage -- purge_haplotigs is
// designed and documented primarily around long-read coverage, so this is
// logged as a caveat rather than presented as an equivalent use case.
//

include { SAMTOOLS_DEPTH                   } from '../../modules/nf-core/samtools/depth/main'
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

    SAMTOOLS_DEPTH(ch_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] })

    PURGEHAPLOTIGS_ESTIMATE_CUTOFFS(SAMTOOLS_DEPTH.out.tsv)
    ch_versions = ch_versions.mix(PURGEHAPLOTIGS_ESTIMATE_CUTOFFS.out.versions)

    PURGEHAPLOTIGS_ESTIMATE_CUTOFFS.out.cutoffs
        .map { meta, cutoffs_json -> [ meta, new groovy.json.JsonSlurper().parse(cutoffs_json) ] }
        .branch { meta, parsed ->
            skip:    parsed.skipped
            proceed: !parsed.skipped
        }
        .set { ch_cutoffs_branch }

    ch_caveats = ch_caveats.mix(
        ch_cutoffs_branch.skip.map { meta, parsed ->
            [ meta, "purge_haplotigs skipped for sample '${meta.id}': ${parsed.reason}" ]
        }
    )
    ch_caveats = ch_caveats.mix(
        ch_cutoffs_branch.proceed.map { meta, parsed ->
            [ meta, "purge_haplotigs cross-check for '${meta.id}' was run on the same coverage used for purge_dups (typically short-read Illumina): purge_haplotigs is designed and documented primarily for long-read coverage, so this should be considered an indicative cross-check, not equivalent to its standard use case." ]
        }
    )

    ch_bam
        .join(ch_fasta)
        .set { ch_for_hist }

    PURGEHAPLOTIGS_HIST(ch_for_hist)
    ch_versions = ch_versions.mix(PURGEHAPLOTIGS_HIST.out.versions)

    PURGEHAPLOTIGS_HIST.out.gencov
        .join(ch_cutoffs_branch.proceed)
        .map { meta, gencov, parsed ->
            [ meta, gencov, parsed.low, parsed.mid, parsed.high ]
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
