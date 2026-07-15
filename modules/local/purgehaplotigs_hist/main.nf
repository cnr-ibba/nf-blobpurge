process PURGEHAPLOTIGS_HIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.3--hdfd78af_0'
        : 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    // purge_haplotigs hist embeds its -d/-depth cutoff (default 200) into
    // both output filenames (e.g. "<bam>.200.gencov",
    // "<bam>.histogram.200.png") -- passed explicitly here so the expected
    // output filenames are documented rather than relying on the tool's
    // default staying 200.
    output:
    tuple val(meta), path("${meta.id}.bam.200.gencov"),          emit: gencov
    tuple val(meta), path("${meta.id}.bam.histogram.200.png"),   emit: histogram_png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ln -s ${bam} ${meta.id}.bam
    ln -s ${bai} ${meta.id}.bam.bai
    samtools faidx ${fasta}

    purge_haplotigs hist -b ${meta.id}.bam -g ${fasta} -t ${task.cpus} -d 200

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_haplotigs: \$( (purge_haplotigs version 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam.200.gencov ${meta.id}.bam.histogram.200.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEHAPLOTIGS_HIST"
    END_VERSIONS
    """
}
