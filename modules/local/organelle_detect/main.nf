process ORGANELLE_DETECT {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    // Classifies contigs as organelle-like from per-contig coverage (+
    // optional reference-match corroboration), and splits the assembly into
    // nuclear/organelle FASTAs -- see bin/detect_organelles.py.
    input:
    tuple val(meta), path(fasta), path(coverage_table), path(reference_paf)

    output:
    tuple val(meta), path("${meta.id}.organelle.fasta"),           emit: organelle_fasta
    tuple val(meta), path("${meta.id}.nuclear.fasta"),             emit: nuclear_fasta
    tuple val(meta), path("${meta.id}.organelle_ids.txt"),         emit: organelle_ids
    tuple val(meta), path("${meta.id}.nuclear_ids.txt"),           emit: nuclear_ids
    tuple val(meta), path("${meta.id}.organelle_report.json"),     emit: report_json
    tuple val(meta), path("${meta.id}.organelle_report.mqc.json"), emit: mqc_json
    tuple val(meta), env('CAVEAT'),                                emit: caveat
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_arg     = params.organelle_reference_fasta ? "--reference-paf ${reference_paf}" : ''
    def max_length_arg    = params.organelle_max_length ? "--max-length ${params.organelle_max_length}" : ''
    def require_match_arg = params.organelle_require_reference_match ? '--require-reference-match' : ''
    """
    detect_organelles.py \\
        --sample-id ${meta.id} \\
        --assembly-fasta ${fasta} \\
        --coverage-table ${coverage_table} \\
        ${reference_arg} \\
        --coverage-multiplier ${params.organelle_coverage_multiplier} \\
        ${max_length_arg} \\
        --min-reference-coverage ${params.organelle_min_reference_coverage} \\
        ${require_match_arg} \\
        --output-json ${meta.id}.organelle_report.json \\
        --output-mqc-json ${meta.id}.organelle_report.mqc.json \\
        --output-caveat ${meta.id}.caveat.txt \\
        --output-organelle-ids ${meta.id}.organelle_ids.txt \\
        --output-nuclear-ids ${meta.id}.nuclear_ids.txt \\
        --output-organelle-fasta ${meta.id}.organelle.fasta \\
        --output-nuclear-fasta ${meta.id}.nuclear.fasta

    CAVEAT=\$(cat ${meta.id}.caveat.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.organelle.fasta ${meta.id}.nuclear.fasta ${meta.id}.organelle_ids.txt ${meta.id}.nuclear_ids.txt
    echo '{"sample_id":"${meta.id}","enabled":true,"n_isolated":0}' > ${meta.id}.organelle_report.json
    echo '{}' > ${meta.id}.organelle_report.mqc.json
    CAVEAT="stub"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "ORGANELLE_DETECT"
    END_VERSIONS
    """
}
