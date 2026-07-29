process HAPLOMERGER2_STAGE_A {
    tag "${meta.id}"
    label 'process_high'

    // HaploMerger2 has no bioconda/conda-forge package and no usable
    // biocontainers image (checked directly against bioconda's package API):
    // it is only distributed as a source tarball, downloaded below at
    // runtime. The container must be supplied via --haplomerger2_container
    // (see docs/usage.md for how to build one, e.g. with Seqera Wave from
    // lastz + UCSC kentUtils bioconda packages); there is no conda profile
    // support for this module.
    container "${params.haplomerger2_container}"

    // Stage A: self-align the raw diploid assembly against itself and break
    // mis-joined scaffolds (hm.batchA1-A3). No coverage input -- unlike
    // purge_dups/purge_haplotigs, HaploMerger2 purges by whole-genome
    // self-alignment alone.
    input:
    tuple val(meta), path(fasta)
    path haplomerger2_path

    output:
    tuple val(meta), path("${meta.id}.haplomerger2_A.fasta.gz"), emit: fasta
    tuple val(meta), path("${meta.id}.haplomerger2_stageA.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def hm2_release   = 'HaploMerger2_20180603'
    def hm2_url       = "https://github.com/mapleforest/HaploMerger2/releases/download/HaploMerger2_20161205/${hm2_release}.tar.gz"
    def use_local_hm2 = haplomerger2_path ? true : false
    def hm2_dir       = use_local_hm2 ? "${haplomerger2_path}" : hm2_release
    """
    if [ "${use_local_hm2}" = "true" ]; then
        HM2_DIR="${hm2_dir}"
    else
        curl -fsSL "${hm2_url}" -o hm2.tar.gz
        tar xzf hm2.tar.gz
        HM2_DIR="${hm2_dir}"
    fi

    # faToNib has no bioconda package; reuse HaploMerger2's own bundled
    # (2010-vintage, statically-built) copy rather than compiling it --
    # everything else (lastz, UCSC kentUtils chain/net tools) is expected on
    # PATH already, from the container built per docs/usage.md.
    export PATH="\$(pwd)/\${HM2_DIR}/chainNet_jksrc20100603_centOS5:\${PATH}"

    # Diploid assemblies with many scaffolds can exceed the default open
    # file-handle limit during self-alignment (HaploMerger2's own docs
    # recommend raising it); best-effort, non-fatal if the runtime forbids it.
    ulimit -n 65535 2>/dev/null || true

    mkdir -p bin project
    cp \${HM2_DIR}/bin/*.pl bin/
    cp \${HM2_DIR}/project_template/*.ctl \${HM2_DIR}/project_template/*.q \${HM2_DIR}/project_template/hm.batchA* project/
    chmod +x project/hm.batchA*

    # identity/threads are hardcoded shell variables inside HaploMerger2's own
    # batch scripts, not CLI flags -- override them in place before running.
    sed -i "s/^identity=.*/identity=${params.haplomerger2_identity}/" project/hm.batchA1.initiation_and_all_lastz
    sed -i "s/^threads=.*/threads=${task.cpus}/"                     project/hm.batchA1.initiation_and_all_lastz
    sed -i "s/^threads=.*/threads=${task.cpus}/"                     project/hm.batchA2.chainNet_and_netToMaf

    zcat -f ${fasta} | gzip -c > project/${meta.id}.fa.gz

    cd project
    ./hm.batchA1.initiation_and_all_lastz ${meta.id}  > ../${meta.id}.haplomerger2_stageA.log 2>&1
    ./hm.batchA2.chainNet_and_netToMaf    ${meta.id} >> ../${meta.id}.haplomerger2_stageA.log 2>&1
    ./hm.batchA3.misjoin_processing       ${meta.id} >> ../${meta.id}.haplomerger2_stageA.log 2>&1
    cd ..

    cp project/${meta.id}_A.fa.gz ${meta.id}.haplomerger2_A.fasta.gz

    # hm.batchA3's own internal pipeline (perl | perl | gzip) has no
    # errexit/pipefail, so a crash partway through (e.g. HM_pathFinder*.pl
    # dying on a too-sparse alignment graph) still leaves this script exiting
    # 0 with an empty/corrupt output -- catch that here instead of silently
    # emitting garbage that would fail confusingly several steps downstream.
    if ! zcat ${meta.id}.haplomerger2_A.fasta.gz 2>/dev/null | grep -q '^>'; then
        echo "ERROR: HaploMerger2 Stage A (hm.batchA1-A3) produced no usable sequence for sample '${meta.id}'." >&2
        echo "This is usually HM_pathFinder_preparation.pl/HM_pathFinder.pl failing on a too-sparse self-alignment graph (too few/short scaffolds), not a configuration error." >&2
        echo "Inspect ${meta.id}.haplomerger2_stageA.log and, in the task work dir, project/*.result/_A3.*.log for the root cause." >&2
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
    printf '>ctg1\\nACGT\\n' | gzip -c > ${meta.id}.haplomerger2_A.fasta.gz
    touch ${meta.id}.haplomerger2_stageA.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "HAPLOMERGER2_STAGE_A"
    END_VERSIONS
    """
}
