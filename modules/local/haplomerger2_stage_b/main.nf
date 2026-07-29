process HAPLOMERGER2_STAGE_B {
    tag "${meta.id}"
    label 'process_high'

    // See HAPLOMERGER2_STAGE_A: no bioconda/conda-forge package, no conda
    // profile support, container supplied via --haplomerger2_container.
    container "${params.haplomerger2_container}"

    // Stage B: self-align the misjoin-corrected assembly from Stage A again,
    // and separate alleles into a reference (kept) and alternative
    // (haplotig) haploid assembly (hm.batchB1-B5).
    input:
    tuple val(meta), path(fasta)
    path haplomerger2_path

    output:
    tuple val(meta), path("${meta.id}.haplomerger2.purged.fasta"),    emit: purged_fasta
    tuple val(meta), path("${meta.id}.haplomerger2.haplotigs.fasta"), emit: removed_fasta
    tuple val(meta), path("${meta.id}.haplomerger2_stageB.log"),      emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def hm2_release   = 'HaploMerger2_20180603'
    def hm2_url       = "https://github.com/mapleforest/HaploMerger2/releases/download/HaploMerger2_20161205/${hm2_release}.tar.gz"
    def use_local_hm2 = haplomerger2_path ? true : false
    def hm2_dir       = use_local_hm2 ? "${haplomerger2_path}" : hm2_release
    def name          = "${meta.id}_A"
    """
    if [ "${use_local_hm2}" = "true" ]; then
        HM2_DIR="${hm2_dir}"
    else
        curl -fsSL "${hm2_url}" -o hm2.tar.gz
        tar xzf hm2.tar.gz
        HM2_DIR="${hm2_dir}"
    fi

    export PATH="\$(pwd)/\${HM2_DIR}/chainNet_jksrc20100603_centOS5:\${PATH}"
    ulimit -n 65535 2>/dev/null || true

    mkdir -p bin project
    cp \${HM2_DIR}/bin/*.pl bin/
    cp \${HM2_DIR}/project_template/*.ctl \${HM2_DIR}/project_template/*.q \${HM2_DIR}/project_template/hm.batchB* project/
    chmod +x project/hm.batchB*

    sed -i "s/^identity=.*/identity=${params.haplomerger2_identity}/" project/hm.batchB1.initiation_and_all_lastz
    sed -i "s/^threads=.*/threads=${task.cpus}/"                     project/hm.batchB1.initiation_and_all_lastz
    sed -i "s/^threads=.*/threads=${task.cpus}/"                     project/hm.batchB4.refine_unpaired_sequences

    zcat -f ${fasta} | gzip -c > project/${name}.fa.gz

    cd project
    ./hm.batchB1.initiation_and_all_lastz            ${name}  > ../${meta.id}.haplomerger2_stageB.log 2>&1
    ./hm.batchB2.chainNet_and_netToMaf               ${name} >> ../${meta.id}.haplomerger2_stageB.log 2>&1
    ./hm.batchB3.haplomerger                         ${name} >> ../${meta.id}.haplomerger2_stageB.log 2>&1
    ./hm.batchB4.refine_unpaired_sequences           ${name} >> ../${meta.id}.haplomerger2_stageB.log 2>&1
    ./hm.batchB5.merge_paired_and_unpaired_sequences ${name} >> ../${meta.id}.haplomerger2_stageB.log 2>&1
    cd ..

    zcat project/${name}_ref.fa.gz > ${meta.id}.haplomerger2.purged.fasta
    zcat project/${name}_alt.fa.gz > ${meta.id}.haplomerger2.haplotigs.fasta

    # hm.batchB3-B5's own internal pipelines have no errexit/pipefail, so a
    # crash partway through (e.g. HM_pathFinder*.pl/XHM_haploMerger.pl dying
    # on a too-sparse alignment graph) still leaves these scripts exiting 0
    # with empty/corrupt output. An empty *reference* assembly is never a
    # legitimate outcome (unlike an empty haplotigs file, which just means no
    # haplotigs were found) -- catch it here instead of failing confusingly
    # several steps downstream (e.g. in BUSCO).
    if ! grep -q '^>' ${meta.id}.haplomerger2.purged.fasta; then
        echo "ERROR: HaploMerger2 Stage B (hm.batchB1-B5) produced no usable reference assembly for sample '${meta.id}'." >&2
        echo "This is usually HM_pathFinder_preparation.pl/HM_pathFinder.pl/XHM_haploMerger.pl failing on a too-sparse self-alignment graph (too few/short scaffolds), not a configuration error." >&2
        echo "Inspect ${meta.id}.haplomerger2_stageB.log and, in the task work dir, project/*.result/_B3.*.log for the root cause." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplomerger2: \$(echo "${hm2_release}" | sed 's/HaploMerger2_//')
        lastz: \$( (lastz --version 2>&1 || true) | grep -m1 -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    echo ">ctg1" > ${meta.id}.haplomerger2.purged.fasta
    echo "ACGT" >> ${meta.id}.haplomerger2.purged.fasta
    touch ${meta.id}.haplomerger2.haplotigs.fasta
    touch ${meta.id}.haplomerger2_stageB.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "HAPLOMERGER2_STAGE_B"
    END_VERSIONS
    """
}
