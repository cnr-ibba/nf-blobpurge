process BUSCO {
    tag "${meta.id}:${stage}:${lineage}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/busco:5.8.3--pyhdfd78af_1'
        : 'quay.io/biocontainers/busco:5.8.3--pyhdfd78af_1'}"

    // Always run with an explicit lineage: --auto-lineage is never used, so
    // results are comparable across stages and reproducible run to run.
    input:
    tuple val(meta), val(stage), path(fasta), val(lineage)
    path busco_lineages_path

    output:
    tuple val(meta), val(stage), val(lineage), path("${meta.id}.${stage}.${lineage}.short_summary.json"), emit: short_summary
    tuple val(meta), val(stage), val(lineage), path("${meta.id}.${stage}.${lineage}_busco"),               emit: full_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}.${stage}.${lineage}"
    def offline_args = busco_lineages_path ? "--offline --download_path ${busco_lineages_path}" : ''
    """
    busco \\
        -i ${fasta} \\
        -o ${prefix} \\
        -l ${lineage} \\
        -m genome \\
        --cpu ${task.cpus} \\
        --out_path . \\
        ${offline_args}

    mv ${prefix}/short_summary.*.json ${prefix}.short_summary.json
    mv ${prefix} ${prefix}_busco

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>/dev/null | sed 's/BUSCO //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${stage}.${lineage}"
    """
    mkdir ${prefix}_busco
    cat <<-END_JSON > ${prefix}.short_summary.json
    {"results": {"Complete percentage": 95.0, "Single copy percentage": 90.0, "Multi copy percentage": 5.0, "Fragmented percentage": 2.0, "Missing percentage": 3.0, "n_markers": 100, "Complete": 95, "Multi copy": 5, "dataset": "${lineage}"}}
    END_JSON
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "BUSCO"
    END_VERSIONS
    """
}
