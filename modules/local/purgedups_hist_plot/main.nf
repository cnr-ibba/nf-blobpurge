process PURGEDUPS_HIST_PLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/matplotlib:3.5.1'
        : 'quay.io/biocontainers/matplotlib:3.5.1'}"

    // The bioconda/biocontainers purge_dups package only ships the compiled
    // C binaries -- purge_dups' own hist_plot.py (vendored into bin/) needs
    // matplotlib, so it runs here rather than in the purge_dups container.
    // Saved as an output (not just the numeric cutoffs) so the calcuts
    // thresholds can be visually sanity-checked against the coverage
    // histogram, as recommended by the purge_dups authors.
    input:
    tuple val(meta), path(stat), path(cutoffs)

    output:
    tuple val(meta), path("${meta.id}.hist.png"), emit: histogram
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    hist_plot.py -c ${cutoffs} ${stat} ${meta.id}.hist.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.hist.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_HIST_PLOT"
    END_VERSIONS
    """
}
