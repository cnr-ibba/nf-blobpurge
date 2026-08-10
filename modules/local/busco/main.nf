process BUSCO {
    tag "${meta.id}:${stage}:${chunk_label}:${lineage}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/busco:5.8.3--pyhdfd78af_1'
        : 'quay.io/biocontainers/busco:5.8.3--pyhdfd78af_1'}"

    // Always run with an explicit lineage: --auto-lineage is never used, so
    // results are comparable across stages and reproducible run to run.
    // Every (meta, stage) fasta is pre-split into span-balanced chunks by
    // SPLIT_ASSEMBLY_FOR_BUSCO upstream (chunk_label is "1of1" when
    // splitting is disabled/not needed) -- chunk_label only disambiguates
    // this task's output filenames; per-chunk results are pooled back into
    // one score per (meta, stage, lineage) by MERGE_BUSCO_CHUNKS.
    input:
    tuple val(meta), val(stage), path(fasta), val(chunk_label), val(lineage)
    path busco_lineages_path

    output:
    tuple val(meta), val(stage), val(lineage), val(chunk_label), path("${meta.id}.${stage}.${chunk_label}.${lineage}.short_summary.json"), emit: short_summary
    tuple val(meta), val(stage), val(lineage), val(chunk_label), path("short_summary.${meta.id}.${stage}.${chunk_label}.${lineage}.txt"),  emit: short_summary_txt
    tuple val(meta), val(stage), val(lineage), val(chunk_label), path("${meta.id}.${stage}.${chunk_label}.${lineage}_busco"),               emit: full_output
    tuple val(meta), val(stage), val(lineage), val(chunk_label), path("${meta.id}.${stage}.${chunk_label}.${lineage}.full_table.tsv"),       emit: full_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}.${stage}.${chunk_label}.${lineage}"
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
    cp ${prefix}/short_summary.*.txt short_summary.${prefix}.txt
    cp ${prefix}/run_${lineage}/full_table.tsv ${prefix}.full_table.tsv
    mv ${prefix} ${prefix}_busco

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>/dev/null | sed 's/BUSCO //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${stage}.${chunk_label}.${lineage}"
    """
    mkdir ${prefix}_busco
    printf '1at100\\tComplete\\tchunk_seq\\t1\\t100\\t+\\t50\\t100\\n' > ${prefix}.full_table.tsv
    cat <<-END_JSON > ${prefix}.short_summary.json
    {"results": {"Complete percentage": 95.0, "Single copy percentage": 90.0, "Multi copy percentage": 5.0, "Fragmented percentage": 2.0, "Missing percentage": 3.0, "n_markers": 100, "Complete": 95, "Multi copy": 5, "dataset": "${lineage}"}}
    END_JSON
    cat <<-END_TXT > short_summary.${prefix}.txt
    # BUSCO version is: 5.8.3
    # The lineage dataset is: ${lineage} (Creation date: stub, number of genomes: 1, number of BUSCOs: 100)

        ***** Results: *****

        95      Complete BUSCOs (C)
        90      Complete and single-copy BUSCOs (S)
        5       Complete and duplicated BUSCOs (D)
        2       Fragmented BUSCOs (F)
        3       Missing BUSCOs (M)
        100     Total BUSCO groups searched
    END_TXT
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "BUSCO"
    END_VERSIONS
    """
}
