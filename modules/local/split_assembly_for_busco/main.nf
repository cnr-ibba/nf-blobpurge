process SPLIT_ASSEMBLY_FOR_BUSCO {
    tag "${meta.id}:${stage}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), val(stage), path(fasta)

    output:
    tuple val(meta), val(stage), path("${meta.id}.${stage}.chunk*.fasta"),      emit: chunks
    tuple val(meta), val(stage), path("${meta.id}.${stage}.split_summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    split_assembly_for_busco.py ${fasta} \\
        --sample-id ${meta.id} \\
        --stage ${stage} \\
        --max-span ${params.busco_chunk_max_span} \\
        --output-prefix ${meta.id}.${stage} \\
        --output-json ${meta.id}.${stage}.split_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${fasta} ${meta.id}.${stage}.chunk1of1.fasta
    echo '{"sample_id":"${meta.id}","stage":"${stage}","chunk_max_span":${params.busco_chunk_max_span},"n_chunks":1,"total_span":4,"total_contigs":1,"chunks":[{"index":1,"n_contigs":1,"span":4}]}' > ${meta.id}.${stage}.split_summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "SPLIT_ASSEMBLY_FOR_BUSCO"
    END_VERSIONS
    """
}
