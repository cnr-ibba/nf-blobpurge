process ORGANELLE_REFERENCE_ALIGN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/minimap2:2.28--h577a1d6_4'
        : 'quay.io/biocontainers/minimap2:2.28--h577a1d6_4'}"

    // Optional corroborating evidence for organelle-contig classification:
    // aligns the filtered assembly (not yet split into nuclear/organelle)
    // against a user-supplied mitochondrial/plastid reference. Skipped
    // (empty PAF) when no reference is given -- coverage alone still drives
    // classification in bin/detect_organelles.py.
    input:
    tuple val(meta), path(fasta), path(reference)

    output:
    tuple val(meta), path("${meta.id}.organelle_ref.paf"), emit: paf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ "${reference.name}" = "NO_FILE" ]; then
        touch ${meta.id}.organelle_ref.paf
    else
        minimap2 -xasm5 -t ${task.cpus} ${reference} ${fasta} > ${meta.id}.organelle_ref.paf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.organelle_ref.paf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "ORGANELLE_REFERENCE_ALIGN"
    END_VERSIONS
    """
}
