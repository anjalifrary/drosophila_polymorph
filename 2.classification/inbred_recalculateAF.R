
library(data.table)
library(SeqArray)

candidate_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.classed.rds")
mel_gds <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/DGRP2.source_BCM-HGSC.dm6.final.reheadered.primaryChr.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.goodSamps.goodSites.gds")
sim_gds <- seqOpen("/scratch/ejy4bu/drosophila/inbred/sampleLevel_filter/dsim3.signor.combined.norm.gatkfilt.snpgap10.snpsOnly.repeatmasked.wmdust.ann.eff.dm6.sorted.goodSamps.goodSites.gds")

### mel:

candidate_dt[, `:=`(
    af_mel_old  = af_mel,
    maf_mel_old = maf_mel,
    n_samps_mel_old = n_samps_mel
)]

seqResetFilter(mel_gds)

seqSetFilterPos(
    mel_gds,
    chr = candidate_dt$chr_dm6,
    pos = candidate_dt$pos_dm6
)

mel_freq_dt <- data.table(
    chr_dm6 = seqGetData(mel_gds, "chromosome"),
    pos_dm6 = seqGetData(mel_gds, "position"),
    af_mel = seqAlleleFreq(mel_gds, ref.allele = 1L)
)

mel_freq_dt[, maf_mel := pmin(af_mel, 1 - af_mel)]

candidate_dt[
    mel_freq_dt,
    on = .(chr_dm6, pos_dm6),
    `:=`(
        af_mel = i.af_mel,
        maf_mel = i.maf_mel
    )
]

candidate_dt[, .(
    n = .N,
    n_af = sum(!is.na(af_mel)),
    n_maf = sum(!is.na(maf_mel)),
    n_old_af = sum(!is.na(af_mel_old))
)]

candidate_dt[
    !is.na(af_mel_old) & !is.na(af_mel),
    .(
        classification,
        chr_dm6,
        pos_dm6,
        af_mel_old,
        af_mel,
        maf_mel_old,
        maf_mel
    )
][1:20]

### sim:

candidate_dt[, `:=`(
    af_sim_old = af_sim,
    maf_sim_old = maf_sim,
    n_samps_sim_old = n_samps_sim
)]

sim_candidates <- candidate_dt[
    !is.na(ref_sim_dm6) & !is.na(alt_sim_dm6)
]

seqResetFilter(sim_gds)

seqSetFilterPos(
    sim_gds,
    chr = sim_candidates$chr_dm6,
    pos = sim_candidates$pos_dm6
)

sim_alleles_dt <- data.table(
    chr_dm6 = seqGetData(sim_gds, "chromosome"),
    pos_dm6  = seqGetData(sim_gds, "position"),
    allele   = seqGetData(sim_gds, "allele")
)
# verify all biallelic alleles retained:
sim_alleles_dt[lengths(strsplit(allele, ",")) > 2] # good!

sim_freq_dt <- data.table(
    chr_dm6 = seqGetData(sim_gds, "chromosome"),
    pos_dm6  = seqGetData(sim_gds, "position"),
    af_sim   = seqAlleleFreq(sim_gds, ref.allele = 1L)
)

sim_freq_dt[, maf_sim := pmin(af_sim, 1 - af_sim)]

candidate_dt[
    sim_freq_dt,
    on = .(chr_dm6, pos_dm6),
    `:=`(
        af_sim = i.af_sim,
        maf_sim = i.maf_sim
    )
]

candidate_dt[, .(
    n = .N,
    n_af = sum(!is.na(af_sim)),
    n_maf = sum(!is.na(maf_sim)),
    n_old_af = sum(!is.na(af_sim_old))
)]

candidate_dt[
    !is.na(af_sim_old) & !is.na(af_sim),
    .(
        classification,
        chr_dm6,
        pos_dm6,
        af_sim_old,
        af_sim,
        maf_sim_old,
        maf_sim
    )
][1:20]

saveRDS(candidate_dt,"/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.candidatesABFGOPXY.classed.MAF.rds")