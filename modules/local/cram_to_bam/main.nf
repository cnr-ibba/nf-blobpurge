process CRAM_TO_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    // Subsets the reads CRAM (aligned upstream to the *original*, unfiltered
    // assembly) down to the contigs that survived BTK_FILTER, and sorts +
    // indexes the result for purge_dups' ngscstat.
    input:
    tuple val(meta), path(original_assembly), path(reads_cram), path(retained_ids)

    output:
    tuple val(meta), path("${meta.id}.coverage.sorted.bam"), path("${meta.id}.coverage.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools faidx ${original_assembly}
    samtools index -@ ${task.cpus} ${reads_cram}

    samtools view -@ ${task.cpus} -b \\
        -T ${original_assembly} \\
        -o ${meta.id}.subset.bam \\
        ${reads_cram} \\
        \$(tr '\\n' ' ' < ${retained_ids})

    samtools sort -@ ${task.cpus} -o ${meta.id}.coverage.sorted.bam ${meta.id}.subset.bam
    samtools index -@ ${task.cpus} ${meta.id}.coverage.sorted.bam
    rm -f ${meta.id}.subset.bam

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
        stub: "CRAM_TO_BAM"
    END_VERSIONS
    """
}
