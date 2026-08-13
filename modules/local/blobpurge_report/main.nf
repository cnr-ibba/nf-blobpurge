process BLOBPURGE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(stats_files), path(span_check), path(busco_comparison_json), path(genomescope_summary), path(organelle_report), val(caveats)

    output:
    tuple val(meta), path("${meta.id}_blobpurge.html"),         emit: html
    tuple val(meta), path("${meta.id}_blobpurge.summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genomescope_arg = genomescope_summary.name != 'NO_FILE' ? "--genomescope-summary ${genomescope_summary}" : ''
    def organelle_arg = organelle_report.name != 'NO_FILE_ORGANELLE' ? "--organelle-report ${organelle_report}" : ''
    def run_ph_arg = params.run_purge_haplotigs ? '--run-purge-haplotigs' : ''
    def caveat_args = caveats.collect { "--caveats '${it.replace("'", "'\\''")}'" }.join(' ')
    """
    generate_report.py \\
        --sample-id ${meta.id} \\
        --stats ${stats_files} \\
        --span-check ${span_check} \\
        --busco-comparison ${busco_comparison_json} \\
        ${genomescope_arg} \\
        ${organelle_arg} \\
        ${run_ph_arg} \\
        ${caveat_args} \\
        --dup-drop-threshold ${params.busco_dup_drop_threshold} \\
        --disagreement-threshold ${params.purge_disagreement_threshold} \\
        --output-html ${meta.id}_blobpurge.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_blobpurge.html
    echo '{}' > ${meta.id}_blobpurge.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "BLOBPURGE_REPORT"
    END_VERSIONS
    """
}
