process FILTER_MIN_CONTIG_LENGTH {
    tag "${meta.id}:${stage}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), val(stage), path(fasta)

    output:
    tuple val(meta), val(stage), path("${meta.id}.${stage}.length_filtered.fasta"), emit: fasta
    tuple val(meta), val(stage), path("${meta.id}.${stage}.contig_length_filter.json"), emit: summary
    tuple val(meta), val(stage), path("${meta.id}.${stage}.contig_length_filter_mqc.json"), emit: mqc_json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_min_contig_length.py ${fasta} \\
        --sample-id ${meta.id} \\
        --stage ${stage} \\
        --min-length ${params.busco_min_contig_length} \\
        --output-fasta ${meta.id}.${stage}.length_filtered.fasta \\
        --output-json ${meta.id}.${stage}.contig_length_filter.json \\
        --output-mqc-json ${meta.id}.${stage}.contig_length_filter_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${fasta} ${meta.id}.${stage}.length_filtered.fasta
    echo '{"sample_id":"${meta.id}","stage":"${stage}","min_contig_length":${params.busco_min_contig_length},"n_contigs_in":1,"n_contigs_out":1,"n_contigs_removed":0,"span_in":4,"span_out":4,"bp_removed":0,"pct_contigs_removed":0.0,"pct_bp_removed":0.0}' > ${meta.id}.${stage}.contig_length_filter.json
    echo '{"id":"contig_length_filter","plot_type":"table","data":{"${meta.id}_${stage}":{"Sample":"${meta.id}","Stage":"${stage}"}}}' > ${meta.id}.${stage}.contig_length_filter_mqc.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "FILTER_MIN_CONTIG_LENGTH"
    END_VERSIONS
    """
}
