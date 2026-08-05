process COMPARE_BUSCO {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(summaries)
    // summaries: BUSCO short_summary*.json files named
    // <sample_id>.<stage>.<lineage>.short_summary.json (see modules/local/busco)

    output:
    tuple val(meta), path("${meta.id}.busco_comparison.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}.busco_comparison.json"), emit: json
    tuple val(meta), path("${meta.id}.busco_comparison_mqc.json"), emit: mqc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    compare_busco.py \\
        --sample-id ${meta.id} \\
        ${summaries} \\
        --dup-drop-threshold ${params.busco_dup_drop_threshold} \\
        --output-tsv ${meta.id}.busco_comparison.tsv \\
        --output-json ${meta.id}.busco_comparison.json \\
        --output-mqc-json ${meta.id}.busco_comparison_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.busco_comparison.tsv
    echo '{"sample_id":"${meta.id}","rows":[],"duplication_drop":{}}' > ${meta.id}.busco_comparison.json
    echo '{"id":"busco_comparison","plot_type":"table","data":{}}' > ${meta.id}.busco_comparison_mqc.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "COMPARE_BUSCO"
    END_VERSIONS
    """
}
