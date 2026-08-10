process ORGANELLE_MERGE_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    // Re-merges organelle contigs (isolated before PURGE_DUPS/PURGE_HAPLOTIGS
    // so they don't distort coverage-cutoff estimation) back into each purge
    // stage's FASTA -- they are real assembly content, not contamination,
    // and stay part of the final assembly. Reuses the python:3.12 image
    // already vendored for BLOBPURGE_REPORT purely as a container with `cat`
    // available; no python is actually invoked (the biocontainers python:3.12
    // image is BusyBox-based, not GNU coreutils -- versions.yml reports the
    // BusyBox version accordingly, not a "coreutils" version).
    input:
    tuple val(meta), path(purged_fasta), path(organelle_fasta)

    output:
    tuple val(meta), path("${prefix}.merged.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    """
    cat ${purged_fasta} ${organelle_fasta} > ${prefix}.merged.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: \$(busybox 2>&1 | head -n1 | sed 's/^BusyBox v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.merged.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "ORGANELLE_MERGE_FASTA"
    END_VERSIONS
    """
}
