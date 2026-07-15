process PURGEDUPS_PURGE_DUPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0'
        : 'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0'}"

    input:
    tuple val(meta), path(cutoffs), path(base_cov), path(self_paf)

    output:
    tuple val(meta), path("${meta.id}.dups.bed"),        emit: bed
    tuple val(meta), path("${meta.id}.purge_dups.log"),  emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-2'
    """
    purge_dups ${args} -T ${cutoffs} -c ${base_cov} ${self_paf} \\
        > ${meta.id}.dups.bed 2> ${meta.id}.purge_dups.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: \$( (purge_dups -h 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.dups.bed ${meta.id}.purge_dups.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_PURGE_DUPS"
    END_VERSIONS
    """
}
