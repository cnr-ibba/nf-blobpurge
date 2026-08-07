process PURGEDUPS_CALCUTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0'
        : 'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0'}"

    input:
    tuple val(meta), path(stat)

    output:
    tuple val(meta), path("${meta.id}.cutoffs"),      emit: cutoffs
    tuple val(meta), path("${meta.id}.calcuts.log"),  emit: log
    tuple val(meta), env('SKIP_REASON'),              emit: skip_reason
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # calcuts needs a genuinely bimodal coverage histogram (a "haploid" peak
    # well separated from the "diploid/collapsed" peak) to derive usable
    # cutoffs. Rather than re-deriving its internal peak-finding statistics
    # ourselves, sanity-check what it actually produced: did it exit cleanly,
    # and did it emit more than one distinct threshold value. A non-zero exit
    # or a degenerate (empty/single-valued) cutoffs file both indicate the
    # coverage distribution wasn't usefully bimodal for this sample -- the
    # caller skips purge_dups for this sample rather than trusting garbage
    # cutoffs or crashing the whole run.
    SKIP_REASON=""
    if ! calcuts ${args} ${stat} > ${meta.id}.cutoffs 2> ${meta.id}.calcuts.log; then
        SKIP_REASON="calcuts exited with an error (see ${meta.id}.calcuts.log); the coverage distribution likely lacks a clean bimodal signal."
    elif [ ! -s ${meta.id}.cutoffs ]; then
        SKIP_REASON="calcuts produced an empty cutoffs file."
    elif [ \$(tr -s ' \\t' '\\n' < ${meta.id}.cutoffs | sort -u | wc -l) -lt 2 ]; then
        SKIP_REASON="calcuts produced degenerate cutoffs (no distinct threshold values); the coverage distribution is not usefully bimodal."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: \$( (purge_dups -h 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cutoffs ${meta.id}.calcuts.log
    SKIP_REASON=""
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_CALCUTS"
    END_VERSIONS
    """
}
