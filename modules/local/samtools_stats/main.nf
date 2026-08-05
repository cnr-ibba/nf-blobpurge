process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_0'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_0'}"

    // QC-only side output of the read_coverage BAM (CRAM subset or fresh
    // bwa-mem2 mapping, whichever path was taken): mapping rate, insert size
    // and per-contig read counts, picked up natively by the MultiQC samtools
    // module. Not consumed by any downstream process.
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.stats"),     emit: stats
    tuple val(meta), path("${meta.id}.flagstat"),  emit: flagstat
    tuple val(meta), path("${meta.id}.idxstats"),  emit: idxstats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${meta.id}.stats
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.id}.flagstat
    samtools idxstats ${bam} > ${meta.id}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.stats ${meta.id}.flagstat ${meta.id}.idxstats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "SAMTOOLS_STATS"
    END_VERSIONS
    """
}
