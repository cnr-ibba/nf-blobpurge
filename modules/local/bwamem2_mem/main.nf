process BWAMEM2_MEM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he70b90d_8'
        : 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'}"

    input:
    tuple val(meta), path(index), path(reads_r1), path(reads_r2)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        \$(ls ${index}/*.0123 | sed 's/\\.0123\$//') \\
        ${reads_r1} ${reads_r2} > ${meta.id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "BWAMEM2_MEM"
    END_VERSIONS
    """
}
