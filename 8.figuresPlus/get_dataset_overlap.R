library(data.table)

# load tables (including non-shared and shared variation, both TSP/conv and different)
poolseq_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/poolseq_dest/classed/dest.mel.sim.PoolSeq.SNAPE.001.50.SynMissense.11_2026.nlp.polyMAF5.shared.bothMelSim.classed.rds")
inbred_dt <- readRDS("/project/berglandlab/anjali/drosophila_polymorphism/classification/inbred/classed/dsim3.signor.DGRP2.source_BCM-HGSC.11_2026.shared.bothMelSim.0.05MAFfilter.classed.rds")


################################################################################
### 1. SNP counts and species overlap within each dataset

dataset_summary <- function(dt, dataset_name) {
  mel <- dt[!is.na(variant.id_mel)]
  sim <- dt[!is.na(variant.id_sim)]

  shared.sameSite <- dt[!is.na(variant.id_mel) & !is.na(variant.id_sim)]

  # Number of genes containing shared variants
  n_genes <- uniqueN(shared.sameSite$gene_id_fbgn, na.rm = TRUE)

  data.table(
    dataset = dataset_name,
    candidates_mel = nrow(mel),
    candidates_sim = nrow(sim),
    shared.sameSite = nrow(shared.sameSite),

    mel_shared_pct = 100 * nrow(shared.sameSite) / nrow(mel),
    sim_shared_pct = 100 * nrow(shared.sameSite) / nrow(sim),

    genes_with_shared = n_genes
  )
}

comp <- rbind(
  dataset_summary(poolseq_dt, "PoolSeq"),
  dataset_summary(inbred_dt, "Inbred")
)

comp

################################################################################
### 2. Candidate overlap between poolseq & inbred

pool_mel <- unique(poolseq_dt[!is.na(variant.id_mel), .(chr, pos)])
inbred_mel <- unique(inbred_dt[!is.na(variant.id_mel), .(chr, pos)])

pool_sim <- unique(poolseq_dt[!is.na(variant.id_sim), .(chr, pos)])
inbred_sim <- unique(inbred_dt[!is.na(variant.id_sim), .(chr, pos)])

mel_overlap <- merge(pool_mel, inbred_mel, by = c("chr", "pos"))

sim_overlap <- merge(pool_sim, inbred_sim, by = c("chr", "pos"))

dataset_overlap <- data.table(
  species = c("D. melanogaster", "D. simulans"),

  poolseq_candidates = c(nrow(pool_mel), nrow(pool_sim)),
  inbred_candidates = c(nrow(inbred_mel), nrow(inbred_sim)),

  overlapping_candidates = c(nrow(mel_overlap), nrow(sim_overlap)),

  pct_poolseq = c(
    100 * nrow(mel_overlap) / nrow(pool_mel),
    100 * nrow(sim_overlap) / nrow(pool_sim)
  ),

  pct_inbred = c(
    100 * nrow(mel_overlap) / nrow(inbred_mel),
    100 * nrow(sim_overlap) / nrow(inbred_sim)
  )
)

dataset_overlap

################################################################################
### 3. Shared polymorphism overlap between PoolSeq and inbred

pool_shared <- unique(poolseq_dt[!is.na(variant.id_mel) & !is.na(variant.id_sim),
    .(chr, pos)])

inbred_shared <- unique(inbred_dt[!is.na(variant.id_mel) & !is.na(variant.id_sim),
    .(chr, pos)])

shared_overlap <- merge(pool_shared, inbred_shared, by = c("chr", "pos"))

shared_summary <- data.table(
  poolseq_shared = nrow(pool_shared),
  inbred_shared = nrow(inbred_shared),
  shared_in_both_datasets = nrow(shared_overlap),

  pct_poolseq_shared_found_in_inbred =
    100 * nrow(shared_overlap) / nrow(pool_shared),

  pct_inbred_shared_found_in_poolseq =
    100 * nrow(shared_overlap) / nrow(inbred_shared)
)
shared_summary