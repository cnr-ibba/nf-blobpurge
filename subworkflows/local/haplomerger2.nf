//
// HAPLOMERGER2: independent cross-check of PURGE_DUPS/PURGE_HAPLOTIGS.
//
// Unlike purge_dups/purge_haplotigs, HaploMerger2 purges by whole-genome
// self-alignment (LASTZ + UCSC kentUtils chain/net) rather than read
// coverage, so it needs no coverage input and forks directly off
// BTK_FILTER's output. Stage A (hm.batchA1-A3) self-aligns the raw diploid
// assembly and breaks mis-joined scaffolds; Stage B (hm.batchB1-B5)
// self-aligns the corrected assembly again and separates alleles into a
// reference (kept) and alternative (haplotig) haploid assembly.
//

include { HAPLOMERGER2_STAGE_A } from '../../modules/local/haplomerger2_stage_a/main'
include { HAPLOMERGER2_STAGE_B } from '../../modules/local/haplomerger2_stage_b/main'

workflow HAPLOMERGER2 {

    take:
    ch_fasta // tuple( meta, filtered_fasta )

    main:
    ch_versions = Channel.empty()
    ch_caveats  = Channel.empty()

    ch_caveats = ch_caveats.mix(
        ch_fasta.map { meta, fasta ->
            [ meta, "HaploMerger2 cross-check per '${meta.id}' si basa esclusivamente su un allineamento self-vs-self (LASTZ) dell'assembly, senza alcuna copertura di lettura: va quindi considerato un cross-check indipendente da purge_dups/purge_haplotigs, non un secondo passaggio sullo stesso segnale, ed e' piu' sensibile alla qualita'/contiguita' dell'assembly e ai parametri di scoring di default (--haplomerger2_identity) che alla profondita' di copertura." ]
        }
    )

    ch_haplomerger2_path = params.haplomerger2_path
        ? Channel.fromPath(params.haplomerger2_path, checkIfExists: true).collect()
        : Channel.value([])

    HAPLOMERGER2_STAGE_A(ch_fasta, ch_haplomerger2_path)
    ch_versions = ch_versions.mix(HAPLOMERGER2_STAGE_A.out.versions)

    HAPLOMERGER2_STAGE_B(HAPLOMERGER2_STAGE_A.out.fasta, ch_haplomerger2_path)
    ch_versions = ch_versions.mix(HAPLOMERGER2_STAGE_B.out.versions)

    emit:
    purged_fasta  = HAPLOMERGER2_STAGE_B.out.purged_fasta
    removed_fasta = HAPLOMERGER2_STAGE_B.out.removed_fasta
    stage_a_log   = HAPLOMERGER2_STAGE_A.out.log
    stage_b_log   = HAPLOMERGER2_STAGE_B.out.log
    caveats       = ch_caveats
    versions      = ch_versions
}
