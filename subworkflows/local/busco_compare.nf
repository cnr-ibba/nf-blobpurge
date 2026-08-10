//
// BUSCO_COMPARE: run BUSCO (never --auto-lineage) on every assembly stage x
// every requested lineage, then build one comparative table per sample.
//

include { BUSCO                     } from '../../modules/local/busco/main'
include { COMPARE_BUSCO             } from '../../modules/local/compare_busco/main'
include { FILTER_MIN_CONTIG_LENGTH  } from '../../modules/local/filter_min_contig_length/main'

workflow BUSCO_COMPARE {

    take:
    ch_stage_fastas // tuple( meta, stage, fasta ) -- one entry per assembly stage to evaluate

    main:
    ch_versions = Channel.empty()
    ch_caveats  = Channel.empty()

    ch_lineages = Channel.fromList(params.busco_lineages.tokenize(','))

    ch_busco_lineages_path = params.busco_lineages_path
        ? Channel.fromPath(params.busco_lineages_path, checkIfExists: true).collect()
        : Channel.value([])

    // BUSCO's own single-threaded post-processing of miniprot's candidate
    // alignments can run out of memory on highly fragmented raw assemblies
    // (large numbers of micro-contigs too short to ever hold a complete gene
    // model, but still searched). Length-filter the raw stage only --
    // filtered/purge_dups/purge_haplotigs have already shed most junk
    // contigs upstream. Disabled entirely via --busco_min_contig_length 0.
    if (params.busco_min_contig_length.toInteger() > 0) {
        ch_stage_fastas
            .branch { meta, stage, fasta ->
                raw:   stage == 'raw'
                other: true
            }
            .set { ch_stage_fastas_branch }

        FILTER_MIN_CONTIG_LENGTH(ch_stage_fastas_branch.raw)
        ch_versions = ch_versions.mix(FILTER_MIN_CONTIG_LENGTH.out.versions)

        FILTER_MIN_CONTIG_LENGTH.out.summary
            .map { meta, _stage, summary_json -> [ meta, new groovy.json.JsonSlurper().parse(summary_json) ] }
            .filter { meta, parsed -> parsed.n_contigs_removed > 0 }
            .map { meta, parsed -> [ meta,
                "BUSCO on the 'raw' assembly stage for '${meta.id}' excludes ${parsed.n_contigs_removed} " +
                "contigs (${parsed.bp_removed} bp, ${String.format('%.2f', parsed.pct_bp_removed)}% of the raw span) " +
                "shorter than --busco_min_contig_length (${parsed.min_contig_length} bp), to avoid BUSCO's " +
                "single-threaded post-processing running out of memory on highly fragmented assemblies. " +
                "ASSEMBLY_STATS and the report's raw-stage span/N50 are still computed from the full, " +
                "unfiltered raw assembly." ]
            }
            .set { ch_caveats }

        ch_busco_input_fastas = ch_stage_fastas_branch.other.mix(FILTER_MIN_CONTIG_LENGTH.out.fasta)
        ch_length_filter_mqc  = FILTER_MIN_CONTIG_LENGTH.out.mqc_json.map { _meta, _stage, f -> f }
    } else {
        ch_busco_input_fastas = ch_stage_fastas
        ch_length_filter_mqc  = Channel.empty()
    }

    ch_busco_input_fastas
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
    length_filter_mqc   = ch_length_filter_mqc
    caveats             = ch_caveats
    versions            = ch_versions
}
