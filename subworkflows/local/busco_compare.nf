//
// BUSCO_COMPARE: run BUSCO (never --auto-lineage) on every assembly stage x
// every requested lineage, then build one comparative table per sample.
//

include { BUSCO         } from '../../modules/local/busco/main'
include { COMPARE_BUSCO } from '../../modules/local/compare_busco/main'

workflow BUSCO_COMPARE {

    take:
    ch_stage_fastas // tuple( meta, stage, fasta ) -- one entry per assembly stage to evaluate

    main:
    ch_versions = Channel.empty()

    ch_lineages = Channel.fromList(params.busco_lineages.tokenize(','))

    ch_busco_lineages_path = params.busco_lineages_path
        ? Channel.fromPath(params.busco_lineages_path, checkIfExists: true).collect()
        : Channel.value([])

    ch_stage_fastas
        .combine(ch_lineages)
        .set { ch_busco_input }

    BUSCO(ch_busco_input, ch_busco_lineages_path)
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    BUSCO.out.short_summary
        .map { meta, _stage, _lineage, summary -> [ meta, summary ] }
        .groupTuple()
        .set { ch_summaries_per_sample }

    COMPARE_BUSCO(ch_summaries_per_sample)
    ch_versions = ch_versions.mix(COMPARE_BUSCO.out.versions)

    emit:
    short_summaries     = BUSCO.out.short_summary
    short_summaries_txt = BUSCO.out.short_summary_txt
    comparison_tsv      = COMPARE_BUSCO.out.tsv
    comparison_json     = COMPARE_BUSCO.out.json
    comparison_mqc      = COMPARE_BUSCO.out.mqc
    versions            = ch_versions
}
