library(data.table)
library(SeqArray)
library(Biostrings)


tigris <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_tigris.snps.NS_SYN.annotate.gds"
steph <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_stephensi.snps.NS_SYN.annotate.gds"
pyrrhus <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_pyrrhus.snps.NS_SYN.annotate.gds"
mitch <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_mitchellii.snps.NS_SYN.annotate.gds"

gds.files <- c(
    tigris = tigris,
    steph = steph,
    pyrrhus = pyrrhus,
    mitch = mitch
)

getEff <- function(genofile) {
    # format of EFF field:
    # EFF, ., String, Predicted effects for this variant.Format: 
    # 'Effect ( Effect_Impact | Functional_Class | Codon_Change | 
    # Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | 
    # Gene_Coding | Transcript_ID | Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )'
    annotated_ids <- seqGetData(genofile, "variant.id") 

    eff_all <- seqGetData(genofile, "annotation/info/EFF")
    id <- rep(annotated_ids, times = eff_all$length)

    data <- eff_all$data
    effect <- sub("\\(.*$", "", data)
    eff_split <- tstrsplit(sub("\\)$", "", sub("^[^(]*\\(", "", data)), "|", fixed=TRUE, fill=NA_character_)

    ### codon info
    refCodon <- sub("/.*", "", eff_split[[3]])
    altCodon <- sub(".*/", "", eff_split[[3]])
    isCodon <- function(codon) {
        return (grepl("/", eff_split[[3]], fixed=TRUE) & nchar(codon) == 3)
    }

    ### amino acid info
    hgvs_p <- ifelse(grepl("^p\\.", eff_split[[4]]), sub("/.*", "", eff_split[[4]]), NA_character_)
    pt_format <- "^p\\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2}|\\*)$"
    isProtein <- function() {
        return (!is.na(hgvs_p) & grepl(pt_format, hgvs_p))
    }
    eff_dt <- data.table(
        id = id, 
        effect = effect,
        impact = eff_split[[1]],
        functional_class = eff_split[[2]],
        gene = eff_split[[6]],
        transcript = eff_split[[9]],
        codon_change = eff_split[[3]],
        ref_codon = ifelse(isCodon(refCodon), refCodon, NA_character_),
        alt_codon = ifelse(isCodon(altCodon), altCodon, NA_character_),
        hgvs_p = hgvs_p,
        ref_aa = ifelse(isProtein(), sub(pt_format, "\\1", hgvs_p), NA_character_),
        alt_aa = ifelse(isProtein(), sub(pt_format, "\\3", hgvs_p), NA_character_),
        aa_pos = ifelse(isProtein(), sub(pt_format, "\\2", hgvs_p), NA_character_),
        allele = eff_split[[11]]

    )
    eff_dt[alt_aa=="*", alt_aa := "Ter"]

    eff_dt
}


getAlleles <- function(genofile, id, ref, altList, tol=0) {
    # tol = tolerance 
    # ref = dt$ref
    # altList = dt$altList
    # tol <- 0
    allele.dt <- data.table(
        id = id,
        ref = ref,
        altList = altList, 
        nAlt = lengths(altList), 
        alt1 = vapply(altList, `[`, character(1), 1),
        alt2 = vapply(altList, `[`, character(1), 2),
        ref0_freq = seqAlleleFreq(genofile, ref.allele=0L),
        alt1_freq = seqAlleleFreq(genofile, ref.allele=1L),
        alt2_freq = seqAlleleFreq(genofile, ref.allele=2L)
    )

    allele.dt[is.na(ref0_freq), ref0_freq := 0]
    allele.dt[is.na(alt1_freq), alt1_freq := 0]
    allele.dt[is.na(alt2_freq), alt2_freq := 0]

    allele.dt[, `:=`(
        ref0_present  = ref0_freq > tol,
        alt1_present = alt1_freq > tol,
        alt2_present = alt2_freq > tol
    )]
    
    # alleles 3+ have no frequency column; only trust the site if 0-2 account for everything (1e-8 is small margin for floating point rounding) 
    allele.dt[, verifiable := ((nAlt <= 2 & nAlt > 0) | (ref0_freq + alt1_freq + alt2_freq) > 1 - 1e-8)]

    allele.dt[, nPoly := rowSums(cbind(ref0_present, alt1_present, alt2_present))]
    allele.dt[, alleles_present := lapply(seq_len(.N), function(i) {
        c(
            if (ref0_present[i]) ref[i],
            if (alt1_present[i]) alt1[i],
            if (alt2_present[i]) alt2[i]
        )
    })] 
    allele.dt[, freq_present := lapply(seq_len(.N), function(i) {
        c(
            if (ref0_present[i]) ref0_freq[i],
            if (alt1_present[i]) alt1_freq[i],
            if (alt2_present[i]) alt2_freq[i]
        )
    })] 

    # biallelic is 2 alleles, is verifiable, and each allele is 1 character [ACGT] only
    allele.dt[, biallelicSNP := (
        nPoly == 2 & 
        verifiable == T & 
        vapply(alleles_present, function(x) all(grepl("^[ACGT]$", x)), logical(1))  
    )]

    allele.dt
}


reconstructCodons <- function(snp.dts, gff_file = "/scratch/ejy4bu/dylan_snakes/snakes/data/Crotalus_pyrrhus_annotation_chr8_SVMP_formatted.gff") {
    gff <- fread(cmd = paste("grep -v '^#'", shQuote(gff_file)), 
        sep="\t", header=F, 
        col.names=c("chr", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")
    )

    cds.dt <- gff[feature == "CDS", .(chr = sub("^chr", "", chr), start, end, strand)]

    snp.overlap <- copy(snp.dts)
    snp.overlap[, `:=`(start = pos, end = pos)]

    setkey(cds.dt, chr, start, end)
    setkey(snp.overlap, chr, start, end)

    snp.cds <- foverlaps(snp.overlap, cds.dt,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type = "within",
        nomatch = NA
    )

    snp.dts[snp.cds, strand := i.strand, on = .(chr, pos = pos)]
    ### many SNPs overlap multiple CDS regions but they are all on the same strand
    snp.cds[, .(
        n_matches = .N,
        n_strands = uniqueN(strand),
        strands = paste(unique(strand), collapse = ",")
    ), by = .(species, chr, pos)][n_strands > 1] # returns 0 

    strand.dt <- unique(snp.cds[, .(chr, pos, strand)])
    snp.dts[strand.dt, strand := i.strand, on = .(chr, pos) ]

    # table(snp.dts$strand, useNA = "always")
    
    snp.dts[, snp_pos_in_codon := regexpr("[A-Z]", sub("/.*", "", codon_change))]

    snp.dts[strand == "+" & snp_pos_in_codon %in% 1:3,
        `:=`(
            codon_genomic_start = pos - (snp_pos_in_codon - 1),
            codon_genomic_end   = pos + (3 - snp_pos_in_codon)
        )]

    snp.dts[strand == "-" & snp_pos_in_codon %in% 1:3,
        `:=`(
            codon_genomic_start = pos - (3 - snp_pos_in_codon),
            codon_genomic_end   = pos + (snp_pos_in_codon - 1)
        )]

    snp.dts[ !is.na(codon_genomic_start), codon_id := paste(chr, codon_genomic_start, codon_genomic_end, sep = ":")]
    # genes [g16607 g16608 g16618 ] are missing from GFF file and therefore missing strandedness & codon_id

    snp.dts
}


snp.dts <- rbindlist(lapply(names(gds.files), function(spp) {
    # gds.fn <- tigris
    # spp <- mitch
    # genofile <- seqOpen(mitch)
    genofile <- seqOpen(gds.files[[spp]])
    seqResetFilter(genofile)

    dt <- data.table(
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"), 
        nAlleles = seqGetData(genofile, "$num_allele"), 
        id = seqGetData(genofile, "variant.id"),
        af = seqGetData(genofile, "annotation/info/AF"), 
        alleles = seqGetData(genofile, "allele"), 
        ref = seqGetData(genofile, "$ref"), 
        altList = strsplit(seqGetData(genofile, "$alt"), ",", fixed=TRUE)
    )

    eff.dt <- getEff(genofile)

    allele.dt <- getAlleles(genofile, dt$id, dt$ref, dt$altList)
    snp.dt <- merge(dt, allele.dt[, !c("ref", "altList"), with=FALSE], by = "id", all.x=TRUE)
    # allele.map <- allele.dt[
    #     biallelicSNP == TRUE,
    #     .(allele = unlist(alleles_present)),
    #     by = id
    # ]
    eff.check <- merge(
        eff.dt,
        allele.dt[, .(id, alleles_present)],
        by = "id",
        allow.cartesian = TRUE
    )

    eff.check[, allele_match := mapply(
        function(a, observed) a %chin% observed,
        allele,
        alleles_present
    )]
    eff.dt <- eff.check[allele_match == TRUE, !c("alleles_present")]


    # snp.dt <- merge(snp.dt, allele.map, by = "id")

    snp.dt <- merge(snp.dt, eff.dt, by = c("id"), all.x=TRUE)
    snp.dt[, species := spp]
    seqClose(genofile)

    snp.dt

    # dt[, is_snp := (nchar(ref) == 1 & nchar(alt) == 1)]
    # dt[, is_indel := (nchar(ref) != nchar(alt))]
}), use.names=T, fill = T)

### filter for biallelic only
snp.dts <- snp.dts[biallelicSNP == T]

### reconstruct codon ids from strand, chr, pos [chr:start:end]
snp.dts <- reconstructCodons(snp.dts)

### update codons with accurate pair of codons
snp.dts[, codons_present := lapply(seq_len(.N), function(i) {
    if (is.na(ref_codon[i]) || is.na(snp_pos_in_codon[i])) {
        return(character(0))
    }

    codon <- strsplit(toupper(ref_codon[i]), "")[[1]]
    alleles <- toupper(unlist(alleles_present[[i]]))

    vapply(alleles, function(a) {
        observed_codon <- codon
        observed_codon[snp_pos_in_codon[i]] <- a
        paste0(observed_codon, collapse = "")
    }, character(1))
})]

### update amino acids with accurate pair of amino acids
snp.dts[, aa_present := lapply(codons_present, function(codons) {
    if (length(codons) == 0) return(character(0))

    as.character(Biostrings::translate(
        DNAStringSet(codons)
    ))
})]

### transcript + AA pos codon id:
snp.dts[, transcriptCodonID := fifelse(
    !is.na(transcript) & !is.na(aa_pos),
    paste(transcript, aa_pos, sep = ":"),
    NA_character_
)]



# ### check codon ids 
# check.dt <- unique(
#     snp.dts[
#         !gene %chin% c("g16607", "g16608", "g16618") &
#         !is.na(codon_id) &
#         !is.na(transcriptCodonID),
#         .(gene, transcript, codon_id, transcriptCodonID)
#     ]
# )
# # Genomic codon -> transcript codon
# check.dt[
#     ,
#     .(n_transcriptCodons = uniqueN(transcriptCodonID)),
#     by = .(gene, transcript, codon_id)
# ][n_transcriptCodons > 1]

# # Transcript codon -> genomic codon
# check.dt[
#     ,
#     .(n_genomicCodons = uniqueN(codon_id)),
#     by = .(gene, transcript, transcriptCodonID)
# ][n_genomicCodons > 1]

###############################################################################################

### CLASSIFY 

cols <- c(
    "chr", "pos", "species", 
    "effect", 
    "snp_pos_in_codon", 
    "aa_pos", 
    "alleles_present", "codons_present", "aa_present", "freq_present",
    "gene", "transcript", "codon_id"
)

codon.long <- snp.dts[
    !is.na(codon_id),
    c("transcriptCodonID", cols),
    with = FALSE
]

codon.sites <- unique(codon.long[, .(codon_id, species, chr, pos, snp_pos_in_codon)])

codon.sites_bySpp <- codon.long[,
    .(transcriptCodonID, snp_pos_in_codon, species, chr, pos, effect, 
    alleles_present, codons_present, aa_present, freq_present, 
    gene, transcript, codon_id )
]

setorder(codon.sites_bySpp, transcriptCodonID, snp_pos_in_codon, species)

### get rows with 1 polymorphic site per species 
codon.sites_bySpp[, n_sites_perCodon_perSpecies := uniqueN(snp_pos_in_codon), by = .(transcriptCodonID, species)]

candidateCodons <- codon.sites_bySpp[n_sites_perCodon_perSpecies == 1]
candidateCodons[, n_species := uniqueN(species), by = .(transcriptCodonID, snp_pos_in_codon)]
candidateCodons <- unique(candidateCodons, by = c("transcriptCodonID", "codon_id", "species", "chr", "pos", "snp_pos_in_codon", "n_species"))

saveRDS(candidateCodons, "/project/berglandlab/anjali/dylan_snakes/snakes/results/snakes.snpDT.candidateCodons.rds")

sameSite <- candidateCodons[n_species>=2, ]

get_pair <- function(x) {
    sort(unique(toupper(unlist(x))))
}

sameSite_summary <- sameSite[
    ,
    {
        codon_pairs <- lapply(codons_present, get_pair)
        aa_pairs <- lapply(aa_present, get_pair)
        names(codon_pairs) <- species
        names(aa_pairs) <- species

        # States present in every species
        shared_codons <- Reduce(intersect, codon_pairs)
        shared_aa <- Reduce(intersect, aa_pairs)

        # All unique states across species
        all_codons <- Reduce(union, codon_pairs)
        all_aa <- Reduce(union, aa_pairs)

        .(
            n_species = uniqueN(species),

            species = list(species),

            codon_pairs = list(codon_pairs),
            aa_pairs = list(aa_pairs),

            total_codons = length(all_codons),
            shared_codons = length(shared_codons),
            Nshared_codons = length(all_codons) - length(shared_codons),

            total_aa = length(all_aa),
            shared_aa = length(shared_aa),
            Nshared_aa = length(all_aa) - length(shared_aa),

            same_codons = length(shared_codons) == length(all_codons),
            same_aa = length(shared_aa) == length(all_aa)
        )
    },
    by = .(transcriptCodonID, codon_id, chr, pos, snp_pos_in_codon)
]

sameSite_summary[same_aa == F, class := "different"]
sameSite_summary[same_aa == T & same_codons == F, class := "convergent"]
sameSite_summary[same_aa ==T & same_codons == T, class := "shared"]

saveRDS(sameSite_summary, "/project/berglandlab/anjali/dylan_snakes/snakes/results/snakes.sameSiteCodons.classes.rds")

candidateCodons <- merge(candidateCodons, sameSite_summary[, .(transcriptCodonID, chr, pos, same_codons, same_aa, class)], by=c("transcriptCodonID", "chr", "pos"), all.x=T)

diffSite_candidates <- candidateCodons[!transcriptCodonID %chin% sameSite$transcriptCodonID]
diffSite_candidates <- diffSite_candidates[, if (uniqueN(species) >= 2  && uniqueN(snp_pos_in_codon) >= 2) .SD, by = transcriptCodonID]
diffSite_summary <- diffSite_candidates[
    ,
    {
        codon_pairs <- lapply(codons_present, get_pair)
        aa_pairs <- lapply(aa_present, get_pair)

        names(codon_pairs) <- species
        names(aa_pairs) <- species

        all_codons <- Reduce(union, codon_pairs)
        all_aa <- Reduce(union, aa_pairs)

        shared_codons <- Reduce(intersect, codon_pairs)
        shared_aa <- Reduce(intersect, aa_pairs)

        .(
            n_species = uniqueN(species),
            species = list(species),
            snp_positions = list(snp_pos_in_codon),

            codon_pairs = list(codon_pairs),
            aa_pairs = list(aa_pairs),

            total_codons = length(all_codons),
            shared_codons = length(shared_codons),
            Nshared_codons = length(all_codons) - length(shared_codons),

            total_aa = length(all_aa),
            shared_aa = length(shared_aa),
            Nshared_aa = length(all_aa) - length(shared_aa),

            same_codons = all(vapply(
                codon_pairs[-1],
                function(x) setequal(x, codon_pairs[[1]]),
                logical(1)
            )),

            same_aa = all(vapply(
                aa_pairs[-1],
                function(x) setequal(x, aa_pairs[[1]]),
                logical(1)
            ))
        )
    },
    by = .(transcriptCodonID, codon_id, gene, transcript)
]