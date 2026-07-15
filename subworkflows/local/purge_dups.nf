//
// PURGE_DUPS: split_fa -> self minimap2 -> calcuts (+ hist_plot.py) -> purge_dups -> get_seqs
//

include { PURGEDUPS_SPLIT_FA      } from '../../modules/local/purgedups_split_fa/main'
include { PURGEDUPS_SELF_MINIMAP2 } from '../../modules/local/purgedups_self_minimap2/main'
include { PURGEDUPS_CALCUTS       } from '../../modules/local/purgedups_calcuts/main'
include { PURGEDUPS_HIST_PLOT     } from '../../modules/local/purgedups_hist_plot/main'
include { PURGEDUPS_PURGE_DUPS    } from '../../modules/local/purgedups_purge_dups/main'
include { PURGEDUPS_GET_SEQS      } from '../../modules/local/purgedups_get_seqs/main'

workflow PURGE_DUPS {

    take:
    ch_fasta // tuple( meta, filtered_fasta )
    ch_stat     // tuple( meta, ngscstat.stat )
    ch_base_cov // tuple( meta, ngscstat.base.cov )

    main:
    ch_versions = Channel.empty()

    PURGEDUPS_SPLIT_FA(ch_fasta)
    ch_versions = ch_versions.mix(PURGEDUPS_SPLIT_FA.out.versions)

    PURGEDUPS_SELF_MINIMAP2(PURGEDUPS_SPLIT_FA.out.fasta)
    ch_versions = ch_versions.mix(PURGEDUPS_SELF_MINIMAP2.out.versions)

    PURGEDUPS_CALCUTS(ch_stat)
    ch_versions = ch_versions.mix(PURGEDUPS_CALCUTS.out.versions)

    ch_stat
        .join(PURGEDUPS_CALCUTS.out.cutoffs)
        .set { ch_for_hist_plot }

    PURGEDUPS_HIST_PLOT(ch_for_hist_plot)
    ch_versions = ch_versions.mix(PURGEDUPS_HIST_PLOT.out.versions)

    PURGEDUPS_CALCUTS.out.cutoffs
        .join(ch_base_cov)
        .join(PURGEDUPS_SELF_MINIMAP2.out.paf)
        .set { ch_for_purge }

    PURGEDUPS_PURGE_DUPS(ch_for_purge)
    ch_versions = ch_versions.mix(PURGEDUPS_PURGE_DUPS.out.versions)

    PURGEDUPS_PURGE_DUPS.out.bed
        .join(ch_fasta)
        .set { ch_for_get_seqs }

    PURGEDUPS_GET_SEQS(ch_for_get_seqs)
    ch_versions = ch_versions.mix(PURGEDUPS_GET_SEQS.out.versions)

    emit:
    purged_fasta  = PURGEDUPS_GET_SEQS.out.purged_fasta
    removed_fasta = PURGEDUPS_GET_SEQS.out.removed_fasta
    cutoffs       = PURGEDUPS_CALCUTS.out.cutoffs
    histogram     = PURGEDUPS_HIST_PLOT.out.histogram
    dups_bed      = PURGEDUPS_PURGE_DUPS.out.bed
    versions      = ch_versions
}
