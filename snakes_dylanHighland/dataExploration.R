library(data.table)
library(SeqArray)

tigris <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_tigris.snps.NS_SYN.annotate.gds"
steph <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_stephensi.snps.NS_SYN.annotate.gds"
pyrrhus <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_pyrrhus.snps.NS_SYN.annotate.gds"
mitch <- "/scratch/ejy4bu/dylan_snakes/snakes/data/gds/crotalus_mitchellii.snps.NS_SYN.annotate.gds"

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
        gt = eff_split[[11]]

    )
    eff_dt[alt_aa=="*", alt_aa := "Ter"]

    eff_dt
}

snp.dt <- foreach(gds.fn=fn)%do%{
    gds.fn <- tigris
    genofile <- seqOpen(gds.fn)
    seqResetFilter(genofile)

    snp.dt <- data.table(
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

    # snp.dt <- merge(eff.dt, snp.dt, by="id", all.x=T)

    getSegregatingAlleles(genofile, ref, altList)

    dt[, is_snp := (nchar(ref) == 1 & nchar(alt) == 1)]
    dt[, is_indel := (nchar(ref) != nchar(alt))]
}


is_snp <- nchar(ref) == 1 & nchar(alt) == 1
    is_indel <- nchar(ref) != nchar(alt)

isSNP <- function(ref, alt) {
    ## TRUE if ref and every alt allele are single nucleotides (A/C/G/T).
    ## Handles multiallelic alts ("A,T"); indels, MNPs, "*", "<DEL>", "." etc. return FALSE.

        ref <- toupper(ref)
        alt <- toupper(alt)
        refOK <- grepl("^[ACGT]$", ref)
        altOK <- grepl("^[ACGT](,[ACGT])*$", alt)
        refOK & altOK
    }


dt <- data.table(
    chr = seqGetData(gds_file, "chromosome"), 
    pos = seqGetData(gds_file, "position"),
    nAlleles = seqGetData(gds_file, "$num_allele"),
    id = seqGetData(gds_file, "variant.id"),
    af = seqGetData(gds_file, "annotation/info/AF") 
    # ,
    # maf = seqGetData(gds_file, "annotation/info/MAF"),
    # n_samps = seqGetData(gds_file, "annotation/info/NS")
)


seqClose(tigris)
seqClose(steph)
seqClose(pyrrhus)
seqClose(mitch)