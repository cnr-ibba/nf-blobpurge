process SAMTOOLS_VIEW_SUBSET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    // Subsets a coordinate-sorted BAM down to a positive list of contig IDs
    // (one per line) -- same samtools view + reindex idiom as CRAM_TO_BAM.
    // An empty ID list must produce an empty BAM, not the unfiltered input:
    // `samtools view` with zero region arguments returns everything, so the
    // empty-list case is handled separately with a header-only BAM.
    input:
    tuple val(meta), path(bam), path(bai), path(ids)

    output:
    tuple val(meta), path("${meta.id}.subset.sorted.bam"), path("${meta.id}.subset.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ -s ${ids} ]; then
        samtools view -@ ${task.cpus} -b \\
            -o ${meta.id}.subset.bam \\
            ${bam} \\
            \$(tr '\\n' ' ' < ${ids})
    else
        samtools view -@ ${task.cpus} -b -H -o ${meta.id}.subset.bam ${bam}
    fi

    samtools sort -@ ${task.cpus} -o ${meta.id}.subset.sorted.bam ${meta.id}.subset.bam
    samtools index ${meta.id}.subset.sorted.bam
    rm -f ${meta.id}.subset.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.subset.sorted.bam ${meta.id}.subset.sorted.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "SAMTOOLS_VIEW_SUBSET"
    END_VERSIONS
    """
}
