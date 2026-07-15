process PURGEDUPS_NGSCSTAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0'
        : 'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0'}"

    // ngscstat is purge_dups' Illumina/short-read coverage-statistics path.
    // Per the purge_dups authors, this path is far less exercised than the
    // PacBio (pbcstat) path -- surfaced here explicitly rather than treated
    // as an equivalent, drop-in replacement. It also requires a name-sorted
    // (not coordinate-sorted) BAM: it silently reports all-zero coverage
    // otherwise, since it pairs up mates assuming they are adjacent.
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.ngscstat.stat"),     emit: stat
    tuple val(meta), path("${meta.id}.ngscstat.base.cov"), emit: base_cov
    tuple val(meta), env(CAVEAT),                          emit: caveat
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ngscstat ${bam}
    mv TX.stat ${meta.id}.ngscstat.stat
    mv TX.base.cov ${meta.id}.ngscstat.base.cov

    CAVEAT="Coverage per purge_dups e' stata calcolata con ngscstat (percorso Illumina/short-read): questo percorso e' meno testato dagli sviluppatori di purge_dups rispetto al percorso PacBio standard (pbcstat) e va considerato con maggiore cautela."
    echo "\${CAVEAT}" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: \$( (purge_dups -h 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.ngscstat.stat ${meta.id}.ngscstat.base.cov
    CAVEAT="stub"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_NGSCSTAT"
    END_VERSIONS
    """
}
