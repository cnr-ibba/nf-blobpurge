//
// BUSCO_COMPARE: run BUSCO (never --auto-lineage) on every assembly stage x
// every requested lineage, then build one comparative table per sample.
//

include { BUSCO                     } from '../../modules/local/busco/main'
include { COMPARE_BUSCO             } from '../../modules/local/compare_busco/main'
include { FILTER_MIN_CONTIG_LENGTH  } from '../../modules/local/filter_min_contig_length/main'
include { SPLIT_ASSEMBLY_FOR_BUSCO  } from '../../modules/local/split_assembly_for_busco/main'
include { MERGE_BUSCO_CHUNKS        } from '../../modules/local/merge_busco_chunks/main'

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
            .set { ch_length_filter_caveats }
        ch_caveats = ch_caveats.mix(ch_length_filter_caveats)

        ch_busco_input_fastas = ch_stage_fastas_branch.other.mix(FILTER_MIN_CONTIG_LENGTH.out.fasta)
        ch_length_filter_mqc  = FILTER_MIN_CONTIG_LENGTH.out.mqc_json.map { _meta, _stage, f -> f }
    } else {
        ch_busco_input_fastas = ch_stage_fastas
        ch_length_filter_mqc  = Channel.empty()
    }

    // Regardless of contig fragmentation, BUSCO's single-threaded
    // post-processing memory scales with total genome content searched
    // against the reference protein set -- large enough to OOM even on a
    // clean, unfragmented assembly (confirmed: the 'filtered' stage, never
    // touched by the length filter above, OOMs too). Split every stage's
    // fasta into span-balanced chunks, run BUSCO independently per chunk,
    // and merge the per-ortholog results back together. A single contig
    // longer than --busco_chunk_max_span still becomes its own oversized
    // chunk (contigs are never split mid-sequence). <=0 disables splitting
    // (always exactly one chunk, handled by the script itself -- no
    // separate code path needed here).
    SPLIT_ASSEMBLY_FOR_BUSCO(ch_busco_input_fastas)
    ch_versions = ch_versions.mix(SPLIT_ASSEMBLY_FOR_BUSCO.out.versions)

    SPLIT_ASSEMBLY_FOR_BUSCO.out.summary
        .map { meta, stage, summary_json -> [ meta, stage, new groovy.json.JsonSlurper().parse(summary_json) ] }
        .filter { meta, stage, parsed -> parsed.n_chunks > 1 }
        .map { meta, stage, parsed -> [ meta,
            "BUSCO on the '${stage}' assembly stage for '${meta.id}' was split into ${parsed.n_chunks} " +
            "span-balanced chunks (target ≤${parsed.chunk_max_span} bp each) and the per-ortholog results " +
            "were merged back together, because BUSCO's single-threaded post-processing does not fit in memory " +
            "on the whole assembly at once. Merging pools each chunk's own Complete/Duplicated/Fragmented calls " +
            "per BUSCO ortholog using the same logic BUSCO applies within a single run, but this has not been " +
            "cross-validated bit-for-bit against an unsplit run on this dataset." ]
        }
        .set { ch_split_caveats }
    ch_caveats = ch_caveats.mix(ch_split_caveats)

    SPLIT_ASSEMBLY_FOR_BUSCO.out.chunks
        .flatMap { meta, stage, chunk_fastas ->
            (chunk_fastas instanceof List ? chunk_fastas : [chunk_fastas]).collect { f ->
                def m = (f.name =~ /\.chunk(\d+)of(\d+)\./)
                [ meta, stage, f, "${m[0][1]}of${m[0][2]}" ]
            }
        }
        .combine(ch_lineages)
        .set { ch_busco_input } // [ meta, stage, fasta, chunk_label, lineage ]

    BUSCO(ch_busco_input, ch_busco_lineages_path)
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    BUSCO.out.full_table
        .map { meta, stage, lineage, _chunk_label, full_table -> [ [meta, stage, lineage], full_table ] }
        .groupTuple()
        .map { key, full_tables -> [ key[0], key[1], key[2], full_tables ] }
        .set { ch_full_tables_grouped } // [ meta, stage, lineage, [full_table.tsv, ...] ]

    MERGE_BUSCO_CHUNKS(ch_full_tables_grouped)
    ch_versions = ch_versions.mix(MERGE_BUSCO_CHUNKS.out.versions)

    MERGE_BUSCO_CHUNKS.out.short_summary
        .map { meta, _stage, _lineage, summary -> [ meta, summary ] }
        .groupTuple()
        .set { ch_summaries_per_sample }

    COMPARE_BUSCO(ch_summaries_per_sample)
    ch_versions = ch_versions.mix(COMPARE_BUSCO.out.versions)

    emit:
    short_summaries     = MERGE_BUSCO_CHUNKS.out.short_summary
    short_summaries_txt = MERGE_BUSCO_CHUNKS.out.short_summary_txt
    comparison_tsv      = COMPARE_BUSCO.out.tsv
    comparison_json     = COMPARE_BUSCO.out.json
    comparison_mqc      = COMPARE_BUSCO.out.mqc
    length_filter_mqc   = ch_length_filter_mqc
    caveats             = ch_caveats
    versions            = ch_versions
}
