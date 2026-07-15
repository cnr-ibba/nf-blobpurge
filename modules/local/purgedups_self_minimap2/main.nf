process PURGEDUPS_SELF_MINIMAP2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/minimap2:2.28--h577a1d6_4'
        : 'quay.io/biocontainers/minimap2:2.28--h577a1d6_4'}"

    input:
    tuple val(meta), path(split_fasta)

    output:
    tuple val(meta), path("${meta.id}.split.self.paf.gz"), emit: paf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 -xasm5 -DP -t ${task.cpus} ${split_fasta} ${split_fasta} \\
        | gzip -c > ${meta.id}.split.self.paf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip -c > ${meta.id}.split.self.paf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_SELF_MINIMAP2"
    END_VERSIONS
    """
}
