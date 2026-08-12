process MERGE_BUSCO_DOWNLOADS {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    path(download_dirs, stageAs: 'lineage_dl_*')
    // download_dirs: one BUSCO_DOWNLOAD.out.download_dir ('busco_downloads')
    // per requested lineage. Every one of them is offline-consumable on its
    // own (each is a complete busco_downloads/{information,lineages,
    // placement_files} tree scoped to that single lineage), so a plain
    // union is enough -- lineages/<lineage>_odbXX subdirs never collide
    // across inputs, and any shared information/placement_files content is
    // expected to be identical since every input comes from the same BUSCO
    // version, so overwriting it while merging is harmless. (Deliberately
    // not using `cp -n`: BusyBox's `cp -n <symlinked-dir>/. dest/` silently
    // copies nothing -- verified against this module's own container.)

    output:
    path "busco_downloads", emit: download_dir
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p busco_downloads
    for lineage_dl in lineage_dl_*; do
        cp -r "\${lineage_dl}"/. busco_downloads/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p busco_downloads
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "MERGE_BUSCO_DOWNLOADS"
    END_VERSIONS
    """
}
