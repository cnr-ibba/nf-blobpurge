process BTK_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/blobtoolkit:4.5.3--pyhdfd78af_0'
        : 'quay.io/biocontainers/blobtoolkit:4.5.3--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(assembly), path(blobdir)

    output:
    tuple val(meta), path("${meta.id}.filtered.fasta"),    emit: fasta
    tuple val(meta), path("retained_ids.txt"),             emit: retained_ids
    tuple val(meta), path("excluded_ids.txt"),             emit: excluded_ids
    tuple val(meta), path("${meta.id}.span_check.json"),   emit: span_check
    tuple val(meta), path("${meta.id}.filter_summary.json"), emit: filter_summary
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix        = task.ext.prefix ?: meta.id
    def exclude_taxa  = params.exclude_taxa
    def taxon_field    = params.taxon_field
    def taxrule_param = params.taxrule ?: ''
    def span_tolerance = params.span_tolerance
    """
    TAXRULE="${taxrule_param}"
    if [ -z "\${TAXRULE}" ]; then
        TAXRULE=\$(detect_taxrule.py ${blobdir}/meta.json)
        echo "Auto-detected --taxrule from BlobDir meta.json: \${TAXRULE}"
    fi

    if ! blobtools filter \\
        --param ${taxon_field}--Keys=${exclude_taxa} \\
        --fasta ${assembly} \\
        --taxrule "\${TAXRULE}" \\
        --output ${prefix}_filtered_btk \\
        --suffix filtered \\
        --summary STDOUT \\
        ${blobdir} > ${prefix}.filter_summary.json
    then
        echo "ERROR: blobtools filter failed for sample '${meta.id}'." >&2
        echo "This is usually caused by --taxon_field ('${taxon_field}') or --taxrule ('\${TAXRULE}')" >&2
        echo "not matching a field actually present in the BlobDir. Inspect" >&2
        echo "${blobdir}/meta.json ('fields' and 'settings.taxrule(s)') and see docs/usage.md." >&2
        exit 1
    fi

    # blobtools filter writes the filtered FASTA next to the input (not
    # inside --output, which only holds the filtered BlobDir field JSONs),
    # inserting ".<suffix>" right after the first "." of the input filename
    # (e.g. "assembly.fasta" -> "assembly.filtered.fasta",
    # "assembly.fasta.gz" -> "assembly.filtered.fasta.gz").
    ASSEMBLY_BASENAME=\$(basename "${assembly}")
    FILTERED_FASTA="\${ASSEMBLY_BASENAME%%.*}.filtered.\${ASSEMBLY_BASENAME#*.}"
    if [ ! -s "\${FILTERED_FASTA}" ]; then
        echo "ERROR: blobtools filter did not produce the expected filtered FASTA file '\${FILTERED_FASTA}'." >&2
        exit 1
    fi
    if [[ "\${FILTERED_FASTA}" == *.gz ]]; then
        zcat "\${FILTERED_FASTA}" > ${prefix}.filtered.fasta
    else
        cp "\${FILTERED_FASTA}" ${prefix}.filtered.fasta
    fi

    check_span.py \\
        --original-fasta ${assembly} \\
        --filtered-fasta ${prefix}.filtered.fasta \\
        --blobdir ${blobdir} \\
        --taxon-field ${taxon_field} \\
        --exclude-taxa "${exclude_taxa}" \\
        --tolerance ${span_tolerance} \\
        --sample-id ${meta.id} \\
        --output-json ${prefix}.span_check.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(blobtools --version 2>&1 | sed 's/^.*blobtoolkit v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    echo ">ctg1" > ${prefix}.filtered.fasta
    echo "ACGT" >> ${prefix}.filtered.fasta
    echo "ctg1" > retained_ids.txt
    echo "ctg2" > excluded_ids.txt
    echo '{"sample_id":"${meta.id}","pass":true,"original_span":8,"filtered_span":4,"excluded_span_from_blobdir":4,"reconstructed_span":8,"relative_diff":0.0,"tolerance":0.001}' > ${meta.id}.span_check.json
    echo '{}' > ${prefix}.filter_summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "BTK_FILTER"
    END_VERSIONS
    """
}
