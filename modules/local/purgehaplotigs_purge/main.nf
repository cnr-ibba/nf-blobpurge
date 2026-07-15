process PURGEHAPLOTIGS_PURGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.3--hdfd78af_0'
        : 'quay.io/biocontainers/purge_haplotigs:1.1.3--hdfd78af_0'}"

    input:
    tuple val(meta), path(fasta), path(coverage_stats)

    output:
    tuple val(meta), path("${meta.id}.curated.fasta"),           emit: purged_fasta
    tuple val(meta), path("${meta.id}.curated.haplotigs.fasta"), emit: removed_fasta
    tuple val(meta), path("${meta.id}.curated.*.log"),           emit: log, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    purge_haplotigs purge ${args} \\
        -g ${fasta} \\
        -c ${coverage_stats} \\
        -o ${meta.id}.curated \\
        -t ${task.cpus}

    # purge_haplotigs exits 0 and writes no output at all when it finds no
    # suspect/artefact contigs to act on ("Nothing left to do, exiting..."):
    # a legitimate outcome (nothing to purge), not a failure. Represent it
    # explicitly as an unchanged assembly rather than crashing on the
    # missing declared outputs.
    if [ ! -s ${meta.id}.curated.fasta ]; then
        echo "purge_haplotigs found nothing to purge for sample '${meta.id}': no suspect/artefact contigs were flagged. Emitting the input assembly unchanged." >&2
        cp ${fasta} ${meta.id}.curated.fasta
        touch ${meta.id}.curated.haplotigs.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_haplotigs: \$( (purge_haplotigs version 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    echo ">ctg1" > ${meta.id}.curated.fasta
    echo "ACGT" >> ${meta.id}.curated.fasta
    touch ${meta.id}.curated.haplotigs.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEHAPLOTIGS_PURGE"
    END_VERSIONS
    """
}
