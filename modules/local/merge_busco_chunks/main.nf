process MERGE_BUSCO_CHUNKS {
    tag "${meta.id}:${stage}:${lineage}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), val(stage), val(lineage), path(full_tables)
    // full_tables: that (sample, stage, lineage)'s per-chunk full_table.tsv
    // files, uniquely named by the BUSCO module's chunk-aware prefix (see
    // modules/local/busco/main.nf).

    output:
    tuple val(meta), val(stage), val(lineage), path("${meta.id}.${stage}.${lineage}.short_summary.json"), emit: short_summary
    tuple val(meta), val(stage), val(lineage), path("short_summary.${meta.id}.${stage}.${lineage}.txt"),  emit: short_summary_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}.${stage}.${lineage}"
    """
    merge_busco_chunks.py \\
        ${full_tables} \\
        --sample-id ${meta.id} \\
        --stage ${stage} \\
        --lineage ${lineage} \\
        --output-json ${prefix}.short_summary.json \\
        --output-txt short_summary.${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${stage}.${lineage}"
    """
    cat <<-END_JSON > ${prefix}.short_summary.json
    {"results": {"Complete percentage": 95.0, "Single copy percentage": 90.0, "Multi copy percentage": 5.0, "Fragmented percentage": 2.0, "Missing percentage": 3.0, "n_markers": 100, "Complete": 95, "Multi copy": 5, "dataset": "${lineage}"}}
    END_JSON
    cat <<-END_TXT > short_summary.${prefix}.txt
    # BUSCO version is: merged from 1 chunk(s)
    # The lineage dataset is: ${lineage} (number of BUSCOs: 100)

    	C:95.0%[S:90.0%,D:5.0%],F:2.0%,M:3.0%,n:100
    	95	Complete BUSCOs (C)
    	90	Complete and single-copy BUSCOs (S)
    	5	Complete and duplicated BUSCOs (D)
    	2	Fragmented BUSCOs (F)
    	3	Missing BUSCOs (M)
    	100	Total BUSCO groups searched
    END_TXT
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "MERGE_BUSCO_CHUNKS"
    END_VERSIONS
    """
}
