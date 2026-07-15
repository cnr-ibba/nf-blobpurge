process SAMTOOLS_DEPTH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.depth.tsv"), emit: depth
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools depth -a ${bam} > ${meta.id}.depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    printf "ctg1\\t1\\t10\\n" > ${meta.id}.depth.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "SAMTOOLS_DEPTH"
    END_VERSIONS
    """
}
