process SAMTOOLS_SORT_NAME {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    // purge_dups' ngscstat pairs up mates by scanning the BAM assuming both
    // reads of a pair are adjacent, and silently reports all-zero coverage
    // on a coordinate-sorted BAM. It needs a name-sorted (or aligner output
    // order) BAM instead; this is kept as a separate output from the
    // coordinate-sorted/indexed BAM used everywhere else (e.g. purge_haplotigs).
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.namesorted.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools sort -n -@ ${task.cpus} -o ${meta.id}.namesorted.bam ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.namesorted.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "SAMTOOLS_SORT_NAME"
    END_VERSIONS
    """
}
