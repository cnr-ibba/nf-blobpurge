process ASSEMBLY_STATS {
    tag "${meta.id}:${stage}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), val(stage), path(fasta)

    output:
    tuple val(meta), val(stage), path("${meta.id}.${stage}.stats.json"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    assembly_stats.py ${fasta} \\
        --sample-id ${meta.id} \\
        --stage ${stage} \\
        --output-json ${meta.id}.${stage}.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","stage":"${stage}","n_contigs":1,"total_span":4,"n50":4,"longest_contig":4}' > ${meta.id}.${stage}.stats.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "ASSEMBLY_STATS"
    END_VERSIONS
    """
}
