process ORGANELLE_EXTRACT_READS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    // Exports the paired-end reads that map to the isolated organelle
    // contigs, for use with external organelle-assembly tools (GetOrganelle,
    // MitoHiFi, oatk, ...) -- this pipeline does not assemble organelles
    // itself. Input must be name-sorted for samtools fastq to pair mates
    // correctly (same requirement as PURGEDUPS_NGSCSTAT).
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.organelle_R1.fastq.gz"), path("${meta.id}.organelle_R2.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools fastq -@ ${task.cpus} \\
        -1 ${meta.id}.organelle_R1.fastq.gz \\
        -2 ${meta.id}.organelle_R2.fastq.gz \\
        -0 /dev/null -s /dev/null -n \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    echo -n "" | gzip -c > ${meta.id}.organelle_R1.fastq.gz
    echo -n "" | gzip -c > ${meta.id}.organelle_R2.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "ORGANELLE_EXTRACT_READS"
    END_VERSIONS
    """
}
