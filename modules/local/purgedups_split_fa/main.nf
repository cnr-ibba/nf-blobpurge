process PURGEDUPS_SPLIT_FA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0'
        : 'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.split.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    split_fa ${fasta} > ${meta.id}.split.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: \$( (purge_dups -h 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.split.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_SPLIT_FA"
    END_VERSIONS
    """
}
