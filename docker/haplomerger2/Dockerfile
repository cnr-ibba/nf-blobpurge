# Container for the HAPLOMERGER2_STAGE_A / HAPLOMERGER2_STAGE_B modules.
#
# No bioconda/conda-forge package or biocontainers image exists for
# HaploMerger2 itself (verified against bioconda's package API); it is only
# distributed as a source tarball, which the modules download and run at
# task runtime. This image only needs to provide HaploMerger2's *external*
# dependencies: lastz plus the UCSC kentUtils tools its core purging path
# (hm.batchA1-A3, hm.batchB1-B5 -- misjoin correction and allele separation;
# not the optional scaffolding/gap-closing stages) actually calls, and perl.
#
# One further kentUtils tool, faToNib, has no bioconda package at all; the
# modules fetch it directly from HaploMerger2's own pinned release tarball
# at runtime instead, so it is intentionally not installed here.
#
# procps-ng provides `ps`, which Nextflow itself (not HaploMerger2) requires
# inside every container to collect task resource-usage metrics -- without
# it, tasks fail immediately with "Command 'ps' ... cannot be found" before
# the process script even runs.
#
# Build and publish with wave-cli (see docs/usage.md), e.g.:
#   wave -f docker/haplomerger2/Dockerfile \
#     --freeze --build-repo docker.io/<user>/nf-blobpurge-haplomerger2 --await

FROM mambaorg/micromamba:2

USER root

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        lastz \
        ucsc-axtchain \
        ucsc-chainantirepeat \
        ucsc-chainfilter \
        ucsc-chainmergesort \
        ucsc-chainnet \
        ucsc-chainprenet \
        ucsc-chainsort \
        ucsc-chainsplit \
        ucsc-chainstitchid \
        ucsc-chainswap \
        ucsc-fafilter \
        ucsc-fasize \
        ucsc-fasplit \
        ucsc-netchainsubset \
        ucsc-netfilter \
        ucsc-netsyntenic \
        ucsc-nettoaxt \
        ucsc-nibfrag \
        perl \
        curl \
        ca-certificates \
        procps-ng \
    && micromamba clean --all --yes

# Set PATH explicitly rather than relying on the base image's entrypoint to
# activate the environment, so the tools are on PATH regardless of how the
# container is invoked (e.g. Nextflow's docker executor with a custom -u).
ENV PATH="/opt/conda/bin:${PATH}"
