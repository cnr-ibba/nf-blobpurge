process PURGEHAPLOTIGS_ESTIMATE_CUTOFFS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    // purge_haplotigs normally requires eyeballing the coverage histogram to
    // pick low/mid/high cutoffs by hand. For a generic, non-interactive
    // pipeline these are instead estimated from the sample's own genome-wide
    // depth distribution (see bin/estimate_purgehaplotigs_cutoffs.py) -- not
    // hardcoded to any organism.
    input:
    tuple val(meta), path(depth_tsv)

    output:
    tuple val(meta), path("${meta.id}.cutoffs.json"),   emit: cutoffs
    tuple val(meta), path("${meta.id}.depth_hist.tsv"), emit: histogram
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    estimate_purgehaplotigs_cutoffs.py \\
        ${depth_tsv} \\
        --sample-id ${meta.id} \\
        --output-json ${meta.id}.cutoffs.json \\
        --histogram-tsv ${meta.id}.depth_hist.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","primary_peak":40,"haploid_peak":20,"low":5,"mid":30,"high":80}' > ${meta.id}.cutoffs.json
    touch ${meta.id}.depth_hist.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEHAPLOTIGS_ESTIMATE_CUTOFFS"
    END_VERSIONS
    """
}
