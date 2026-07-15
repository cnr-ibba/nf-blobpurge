process BWAMEM2_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he70b90d_8'
        : 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2_index"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir bwamem2_index
    bwa-mem2 index -p bwamem2_index/${fasta.baseName} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n1)
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2_index
    touch bwamem2_index/${fasta.baseName}.0123
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "BWAMEM2_INDEX"
    END_VERSIONS
    """
}
