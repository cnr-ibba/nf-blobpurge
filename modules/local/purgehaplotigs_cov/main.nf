process PURGEHAPLOTIGS_COV {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.3--hdfd78af_0'
        : 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'}"

    input:
    tuple val(meta), path(gencov), val(low), val(mid), val(high)

    output:
    tuple val(meta), path("${meta.id}.coverage_stats.csv"), emit: coverage_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    purge_haplotigs cov ${args} \\
        -i ${gencov} \\
        -l ${low} -m ${mid} -h ${high} \\
        -o ${meta.id}.coverage_stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_haplotigs: \$( (purge_haplotigs version 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_stats.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEHAPLOTIGS_COV"
    END_VERSIONS
    """
}
