process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.coverage.sorted.bam"), path("${meta.id}.coverage.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.coverage.sorted.bam ${alignment}
    samtools index -@ ${task.cpus} ${meta.id}.coverage.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage.sorted.bam ${meta.id}.coverage.sorted.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "SAMTOOLS_SORT_INDEX"
    END_VERSIONS
    """
}
