process CRAM_TO_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    // Subsets the reads CRAM (aligned upstream to the *original*, unfiltered
    // assembly) down to a positive list of contig IDs, and sorts + indexes
    // the result. An empty retained_ids file is a deliberate "no filtering"
    // sentinel: `samtools view` with zero region arguments returns every
    // contig, so this same module also serves as a raw-assembly-wide BAM
    // producer for ORGANELLE_ISOLATE (fed assets/NO_FILE as retained_ids) --
    // see SAMTOOLS_VIEW_SUBSET for the same idiom used the other way around.
    input:
    tuple val(meta), path(original_assembly), path(reads_cram), path(retained_ids)

    output:
    tuple val(meta), path("*.coverage.sorted.bam"), path("*.coverage.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    samtools faidx ${original_assembly}
    samtools index -@ ${task.cpus} ${reads_cram}

    samtools view -@ ${task.cpus} -b \\
        -T ${original_assembly} \\
        -o ${prefix}.subset.bam \\
        ${reads_cram} \\
        \$(tr '\\n' ' ' < ${retained_ids})

    samtools sort -@ ${task.cpus} -o ${prefix}.coverage.sorted.bam ${prefix}.subset.bam
    samtools index -@ ${task.cpus} ${prefix}.coverage.sorted.bam
    rm -f ${prefix}.subset.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.coverage.sorted.bam ${prefix}.coverage.sorted.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "CRAM_TO_BAM"
    END_VERSIONS
    """
}
