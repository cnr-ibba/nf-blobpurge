process PURGEDUPS_GET_SEQS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0'
        : 'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0'}"

    input:
    tuple val(meta), path(dups_bed), path(filtered_fasta)

    output:
    tuple val(meta), path("${meta.id}.purged.fasta"),          emit: purged_fasta
    tuple val(meta), path("${meta.id}.purged_haplotigs.fasta"), emit: removed_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // No -e: this is a primary (not haplotype-merged) assembly, so partial
    // overlaps are allowed to be trimmed rather than only whole contigs removed.
    def args = task.ext.args ?: ''
    """
    get_seqs ${args} ${dups_bed} ${filtered_fasta}
    mv purged.fa ${meta.id}.purged.fasta
    mv hap.fa ${meta.id}.purged_haplotigs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: \$( (purge_dups -h 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    echo ">ctg1" > ${meta.id}.purged.fasta
    echo "ACGT" >> ${meta.id}.purged.fasta
    touch ${meta.id}.purged_haplotigs.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "PURGEDUPS_GET_SEQS"
    END_VERSIONS
    """
}
