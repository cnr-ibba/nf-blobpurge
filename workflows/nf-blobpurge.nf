/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BTK_FILTER              } from '../modules/local/btk_filter/main'
include { ASSEMBLY_STATS          } from '../modules/local/assembly_stats/main'
include { BLOBPURGE_REPORT        } from '../modules/local/blobpurge_report/main'
include { READ_COVERAGE           } from '../subworkflows/local/read_coverage'
include { ORGANELLE_ISOLATE       } from '../subworkflows/local/organelle_isolate'
include { ORGANELLE_MERGE_FASTA as ORGANELLE_MERGE_FASTA_FILTERED       } from '../modules/local/organelle_merge_fasta/main'
include { ORGANELLE_MERGE_FASTA as ORGANELLE_MERGE_FASTA_PURGEDUPS      } from '../modules/local/organelle_merge_fasta/main'
include { ORGANELLE_MERGE_FASTA as ORGANELLE_MERGE_FASTA_PURGEHAPLOTIGS } from '../modules/local/organelle_merge_fasta/main'
include { PURGE_DUPS              } from '../subworkflows/local/purge_dups'
include { PURGE_HAPLOTIGS         } from '../subworkflows/local/purge_haplotigs'
include { BUSCO_COMPARE           } from '../subworkflows/local/busco_compare'
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_nf-blobpurge_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BLOBPURGE {

    take:
    ch_samplesheet // channel: [ meta, assembly, blobdir, reads_cram ]
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions       = channel.empty()
    def ch_multiqc_files  = channel.empty()
    def ch_stage_fastas   = channel.empty() // [ meta, stage, fasta ]
    def ch_caveats        = channel.empty() // [ meta, caveat_string ]

    ch_stage_fastas = ch_stage_fastas.mix(
        ch_samplesheet.map { meta, assembly, blobdir, reads_cram -> [ meta, 'raw', assembly ] }
    )

    //
    // STEP 0.5: organelle contig isolation (optional), BEFORE contamination
    // filtering. Mitochondria are frequently taxonomically misclassified as
    // bacterial contamination, so classifying organelle-like contigs (high
    // copy number -> extreme coverage) off the RAW assembly + its
    // already-aligned CRAM -- and physically splitting them out before
    // BTK_FILTER's `--fasta` input is built -- means blobtools' taxonomic
    // filter can never see, and so can never discard, an organelle contig.
    // This also keeps them out of purge_dups'/purge_haplotigs' coverage-
    // cutoff estimation downstream, with no separate accounting needed there
    // (BTK_FILTER's own output is organelle-free by construction once its
    // input is). Off by default -- when off, every channel below is a plain
    // alias of today's channels, so the DAG and pipeline behaviour are
    // unchanged.
    //
    def ch_fasta_for_filter = ch_samplesheet.map { meta, assembly, blobdir, reads_cram -> [ meta, assembly ] }
    def ch_organelle_ids    = ch_samplesheet
        .map { meta, assembly, blobdir, reads_cram -> meta }
        .combine(Channel.fromPath("${projectDir}/assets/NO_FILE"))
    def ch_organelle_fasta  = channel.empty() // [ meta, organelle_fasta ], only populated when the feature is on

    if (params.run_organelle_isolation) {
        def ch_organelle_reference = params.organelle_reference_fasta
            ? Channel.fromPath(params.organelle_reference_fasta, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_FILE")

        ch_samplesheet
            .map { meta, assembly, blobdir, reads_cram -> meta }
            .combine(ch_organelle_reference)
            .set { ch_organelle_reference_per_sample }

        ORGANELLE_ISOLATE(
            ch_samplesheet.map { meta, assembly, blobdir, reads_cram -> [ meta, assembly, reads_cram ] },
            ch_organelle_reference_per_sample
        )
        ch_versions = ch_versions.mix(ORGANELLE_ISOLATE.out.versions)
        ch_caveats  = ch_caveats.mix(ORGANELLE_ISOLATE.out.caveats)
        ch_multiqc_files = ch_multiqc_files.mix(ORGANELLE_ISOLATE.out.multiqc_files)

        ch_fasta_for_filter = ORGANELLE_ISOLATE.out.nuclear_fasta
        ch_organelle_ids    = ORGANELLE_ISOLATE.out.organelle_ids
        ch_organelle_fasta  = ORGANELLE_ISOLATE.out.organelle_fasta
    }

    //
    // STEP 1: contamination filtering, cross-checked against the BlobDir
    //
    BTK_FILTER(
        ch_fasta_for_filter
            .join(ch_samplesheet.map { meta, assembly, blobdir, reads_cram -> [ meta, blobdir ] })
            .join(ch_organelle_ids)
    )
    ch_versions = ch_versions.mix(BTK_FILTER.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BTK_FILTER.out.span_check_mqc.map { _meta, f -> f })

    def ch_filtered_stage_fasta = BTK_FILTER.out.fasta
    if (params.run_organelle_isolation) {
        ORGANELLE_MERGE_FASTA_FILTERED(BTK_FILTER.out.fasta.join(ch_organelle_fasta))
        ch_versions = ch_versions.mix(ORGANELLE_MERGE_FASTA_FILTERED.out.versions)
        ch_filtered_stage_fasta = ORGANELLE_MERGE_FASTA_FILTERED.out.fasta
    }
    ch_stage_fastas = ch_stage_fastas.mix(ch_filtered_stage_fasta.map { meta, fasta -> [ meta, 'filtered', fasta ] })

    //
    // STEP 2: read coverage for purge_dups (CRAM subset onto BTK_FILTER-retained contigs)
    //
    ch_samplesheet
        .map { meta, assembly, blobdir, reads_cram -> [ meta, assembly, reads_cram ] }
        .join(BTK_FILTER.out.retained_ids)
        .set { ch_for_coverage }

    READ_COVERAGE(ch_for_coverage)
    ch_versions = ch_versions.mix(READ_COVERAGE.out.versions)
    ch_caveats  = ch_caveats.mix(READ_COVERAGE.out.caveats)
    ch_multiqc_files = ch_multiqc_files.mix(READ_COVERAGE.out.multiqc_files)

    //
    // STEP 3: purge_dups (always run). BTK_FILTER.out.fasta/READ_COVERAGE's
    // outputs are already organelle-free by construction whether or not
    // organelle isolation is on -- no ch_*_for_purge indirection needed.
    //
    PURGE_DUPS(
        BTK_FILTER.out.fasta,
        READ_COVERAGE.out.stat,
        READ_COVERAGE.out.base_cov
    )
    ch_versions = ch_versions.mix(PURGE_DUPS.out.versions)
    ch_caveats  = ch_caveats.mix(PURGE_DUPS.out.caveats)

    def ch_purge_dups_stage_fasta = PURGE_DUPS.out.purged_fasta
    if (params.run_organelle_isolation) {
        ORGANELLE_MERGE_FASTA_PURGEDUPS(PURGE_DUPS.out.purged_fasta.join(ch_organelle_fasta))
        ch_versions = ch_versions.mix(ORGANELLE_MERGE_FASTA_PURGEDUPS.out.versions)
        ch_purge_dups_stage_fasta = ORGANELLE_MERGE_FASTA_PURGEDUPS.out.fasta
    }
    ch_stage_fastas = ch_stage_fastas.mix(ch_purge_dups_stage_fasta.map { meta, fasta -> [ meta, 'purge_dups', fasta ] })

    //
    // STEP 4: purge_haplotigs cross-check (optional, independent, in parallel to purge_dups)
    //
    if (params.run_purge_haplotigs) {
        PURGE_HAPLOTIGS(
            BTK_FILTER.out.fasta,
            READ_COVERAGE.out.bam
        )
        ch_versions = ch_versions.mix(PURGE_HAPLOTIGS.out.versions)
        ch_caveats  = ch_caveats.mix(PURGE_HAPLOTIGS.out.caveats)

        def ch_purge_haplotigs_stage_fasta = PURGE_HAPLOTIGS.out.purged_fasta
        if (params.run_organelle_isolation) {
            ORGANELLE_MERGE_FASTA_PURGEHAPLOTIGS(PURGE_HAPLOTIGS.out.purged_fasta.join(ch_organelle_fasta))
            ch_versions = ch_versions.mix(ORGANELLE_MERGE_FASTA_PURGEHAPLOTIGS.out.versions)
            ch_purge_haplotigs_stage_fasta = ORGANELLE_MERGE_FASTA_PURGEHAPLOTIGS.out.fasta
        }
        ch_stage_fastas = ch_stage_fastas.mix(ch_purge_haplotigs_stage_fasta.map { meta, fasta -> [ meta, 'purge_haplotigs', fasta ] })
    }

    //
    // STEP 5: traceable span/N50 per stage, and comparative BUSCO per stage x lineage
    //
    ASSEMBLY_STATS(ch_stage_fastas)
    ch_versions = ch_versions.mix(ASSEMBLY_STATS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_STATS.out.mqc_json.map { _meta, _stage, f -> f })

    ASSEMBLY_STATS.out.stats
        .map { meta, _stage, stats -> [ meta, stats ] }
        .groupTuple()
        .set { ch_stats_grouped }

    BUSCO_COMPARE(ch_stage_fastas)
    ch_versions = ch_versions.mix(BUSCO_COMPARE.out.versions)
    ch_caveats  = ch_caveats.mix(BUSCO_COMPARE.out.caveats)
    ch_multiqc_files = ch_multiqc_files.mix(
        BUSCO_COMPARE.out.short_summaries_txt.map { _meta, _stage, _lineage, summary -> summary }
    )
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_COMPARE.out.comparison_mqc.map { _meta, f -> f })
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_COMPARE.out.length_filter_mqc)

    //
    // STEP 6: per-sample report
    //
    def ch_genomescope = params.genomescope_summary
        ? Channel.fromPath(params.genomescope_summary, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_FILE")

    ch_samplesheet
        .map { meta, assembly, blobdir, reads_cram -> meta }
        .combine(ch_genomescope)
        .set { ch_genomescope_per_sample }

    ch_caveats
        .groupTuple()
        .map { meta, caveats -> [ meta, caveats.sort() ] }
        .set { ch_caveats_grouped }

    // Distinct sentinel from ch_genomescope's assets/NO_FILE: BLOBPURGE_REPORT
    // takes both in the same input tuple, and Nextflow refuses to stage two
    // different input files under the identical literal name "NO_FILE" into
    // one task's work dir.
    def ch_organelle_report = params.run_organelle_isolation
        ? ORGANELLE_ISOLATE.out.report_json
        : ch_samplesheet
            .map { meta, assembly, blobdir, reads_cram -> meta }
            .combine(Channel.fromPath("${projectDir}/assets/NO_FILE_ORGANELLE"))

    ch_stats_grouped
        .join(BTK_FILTER.out.span_check)
        .join(BUSCO_COMPARE.out.comparison_json)
        .join(ch_genomescope_per_sample)
        .join(ch_organelle_report)
        .join(ch_caveats_grouped, remainder: true)
        .map { meta, stats, span_check, busco_json, genomescope, organelle_report, caveats ->
            // ch_caveats is mixed from several independent, parallel subworkflows
            // (READ_COVERAGE, ORGANELLE_ISOLATE, PURGE_DUPS, PURGE_HAPLOTIGS,
            // BUSCO_COMPARE); their relative completion order -- and therefore
            // groupTuple()'s arrival order -- is not deterministic between runs.
            // Sorting here keeps the report's caveat list (and hence its md5)
            // stable regardless of run-to-run scheduling.
            [ meta, stats, span_check, busco_json, genomescope, organelle_report, (caveats ?: []).sort() ]
        }
        .set { ch_for_report }

    BLOBPURGE_REPORT(ch_for_report)
    ch_versions = ch_versions.mix(BLOBPURGE_REPORT.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'nf-blobpurge_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'nf-blobpurge'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
