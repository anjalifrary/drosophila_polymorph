
### functions
    isSNP <- function(ref, alt) {
        ## TRUE if ref and every alt allele are single nucleotides (A/C/G/T).
        ## Handles multiallelic alts ("A,T"); indels, MNPs, "*", "<DEL>", "." etc. return FALSE.

        ref <- toupper(ref)
        alt <- toupper(alt)
        refOK <- grepl("^[ACGT]$", ref)
        altOK <- grepl("^[ACGT](,[ACGT])*$", alt)
        refOK & altOK
    }

    getBiallelic <- function(ref, alt, ref_0_freq, alt_1_freq, alt_2_freq, tol=0) {
        ## Given ref, alt (comma-separated) and the frequencies of alleles 0, 1, 2,
        ## identify the single alt allele that is actually segregating.
        ## Returns a data.table with ref, alt_used, ref_freq, alt_freq.
        ##   - polymorphic (ref + one alt > tol): observed frequencies
        ##   - fixed for ref:  ref_freq = 1, alt_freq = 0, alt_used = first listed alt
        ##   - fixed for an alt: ref_freq = 0, alt_freq = 1, alt_used = that alt
        ## Rows get NA for alt_used/ref_freq/alt_freq when:
        ##   - all frequencies NaN (no called genotypes)
        ##   - ref + both alts polymorphic, or two alts polymorphic but ref absent
        ##   - alleles 3+ exist and may carry frequency (freqs of alleles 0-2 don't sum to 1)

        f <- cbind(ref_0_freq, alt_1_freq, alt_2_freq)
        f[is.na(f)] <- 0                       # NaN alt_2 at biallelic sites -> allele absent

        poly  <- f > tol
        nPoly <- rowSums(poly)

        altList <- strsplit(alt, ",", fixed=TRUE)
        nAlt    <- lengths(altList)
        alt1    <- vapply(altList, `[`, character(1), 1)
        alt2    <- vapply(altList, `[`, character(1), 2)

        ## alleles 3+ have no frequency column; only trust the site if 0-2 account for everything
        verifiable <- nAlt <= 2 | rowSums(f) > 1 - 1e-8

        seg      <- nPoly == 2 & poly[,1] & verifiable          # ref + one alt segregating
        fixedRef <- nPoly == 1 & poly[,1] & verifiable
        fixedAlt <- nPoly == 1 & !poly[,1] & verifiable

        useAlt2 <- poly[,3] & !poly[,2]                          # alt 2 is the one present
        out <- data.table(ref      = ref,
                          alt_used = ifelse(useAlt2, alt2, alt1),
                          ref_freq = f[,1],
                          alt_freq = ifelse(useAlt2, f[,3], f[,2]))
        out[fixedRef, c("ref_freq", "alt_freq") := list(1, 0)]
        out[fixedAlt, c("ref_freq", "alt_freq") := list(0, 1)]
        out[!(seg | fixedRef | fixedAlt),
            c("alt_used", "ref_freq", "alt_freq") := list(NA_character_, NA_real_, NA_real_)]

        out
    }

    getEFF <- function(genofile) {
        ## Parse the snpEff EFF annotation into one row per effect (a variant can have several:
        ## one per alt allele and/or transcript).
        ## EFF format: Effect(Impact|Class|Codon_Change|AA_Change/cDNA|AA_len|Gene|Biotype|Coding|Transcript|Exon|Allele[|Warnings])
        ## Returns: id (variant.id), eff_allele (genomic alt allele the effect refers to), effect, impact,
        ##          gene, transcript, ref_codon, alt_codon, ref_aa, alt_aa, aa_pos, hgvs_p.
        ## Codons are NA for non-coding effects (intron, intragenic, up/downstream).
        ## ref_aa/alt_aa/aa_pos are filled for single-codon changes (synonymous, missense, stop_gained,
        ## e.g. p.Thr33Ser, p.Gln110*); indels/frameshifts keep only hgvs_p (e.g. p.Gly37dup, p.Asp588fs).

        x  <- seqGetData(genofile, "annotation/info/EFF")
        id <- rep(seqGetData(genofile, "variant.id"), x$length)
        v  <- x$data

        effect <- sub("\\(.*", "", v)
        fields <- tstrsplit(sub("\\)$", "", sub("^[^(]*\\(", "", v)), "|", fixed=TRUE, fill=NA_character_)

        codon  <- fields[[3]]
        isCod  <- grepl("/", codon, fixed=TRUE)
        hgvs_p <- ifelse(grepl("^p\\.", fields[[4]]), sub("/.*", "", fields[[4]]), NA_character_)

        aaPat  <- "^p\\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2}|\\*)$"
        isPt   <- !is.na(hgvs_p) & grepl(aaPat, hgvs_p)

        out <- data.table(id         = id,
                          eff_allele = fields[[11]],
                          effect     = effect,
                          impact     = fields[[1]],
                          gene       = fields[[6]],
                          transcript = fields[[9]],
                          ref_codon  = ifelse(isCod, sub("/.*", "", codon), NA_character_),
                          alt_codon  = ifelse(isCod, sub(".*/", "", codon), NA_character_),
                          ref_aa     = ifelse(isPt, sub(aaPat, "\\1", hgvs_p), NA_character_),
                          alt_aa     = ifelse(isPt, sub(aaPat, "\\3", hgvs_p), NA_character_),
                          aa_pos     = replace(rep(NA_integer_, length(v)), isPt,
                                               as.integer(sub(aaPat, "\\2", hgvs_p[isPt]))),
                          hgvs_p     = hgvs_p)
        out[alt_aa == "*", alt_aa := "Ter"]   # stop codon, to match 3-letter style
        out
    }

    classifyPoly <- function(snp.dt) {
        ## Pairwise comparison of amino-acid polymorphisms found in >= 2 species.
        ## Uses sites polymorphic within a species (0 < ref_freq < 1) with a single-codon change
        ## (synonymous, missense, stop_gained; i.e. aa_pos not NA).
        ## Species are compared at the same amino-acid position (transcript + aa_pos), so
        ## changes at different nucleotide positions of the same codon are compared too.
        ##   shared     : same genomic position, same ref/alt nucleotides, same codon change
        ##                (hence same amino acids and same NS/S class)
        ##   convergent : different nucleotides / codon change, but same pair of amino acids
        ##                (NS: e.g. Thr<->Ile reached two ways; S: e.g. Leu<->Leu via different codons)
        ##   different  : same amino-acid position, different amino-acid changes
        ## Returns one row per species pair per polymorphism pair.

        p <- snp.dt[ref_freq > 0 & ref_freq < 1 & !is.na(aa_pos),
                    .(spp, gene, transcript, aa_pos, chr, pos, ref_used, alt_used,
                      ref_codon = toupper(ref_codon), alt_codon = toupper(alt_codon),
                      ref_aa, alt_aa, type = ifelse(ref_aa == alt_aa, "S", "NS"))]

        pr <- merge(p, p, by=c("gene", "transcript", "aa_pos"), suffixes=c("_1", "_2"),
                    allow.cartesian=TRUE)[spp_1 < spp_2]

        pr[, sameNt    := chr_1 == chr_2 & pos_1 == pos_2 & ref_used_1 == ref_used_2 & alt_used_1 == alt_used_2]
        pr[, sameCodon := ref_codon_1 == ref_codon_2 & alt_codon_1 == alt_codon_2]
        ## amino-acid pair compared unordered, in case ref/alt polarity differs between species
        pr[, sameAA    := pmin(ref_aa_1, alt_aa_1) == pmin(ref_aa_2, alt_aa_2) &
                          pmax(ref_aa_1, alt_aa_1) == pmax(ref_aa_2, alt_aa_2)]

        pr[, class := fifelse(sameNt & sameCodon, "shared",
                      fifelse(sameAA & !sameCodon & !sameNt, "convergent", "different"))]
        pr[, type := fifelse(type_1 == type_2, type_1, "mixed")]

        ## When a codon has >1 polymorphism, every SNP in species 1 is crossed with every SNP in
        ## species 2. If a SNP already has an exact (shared) match in the other species, drop its
        ## cross-comparisons with the other SNPs in that codon - they are not meaningful.
        pr[, shared_1 := any(class == "shared"), .(transcript, aa_pos, spp_1, spp_2, pos_1)]
        pr[, shared_2 := any(class == "shared"), .(transcript, aa_pos, spp_1, spp_2, pos_2)]
        pr <- pr[class == "shared" | !(shared_1 | shared_2)]
        pr[, c("shared_1", "shared_2") := NULL]

        setcolorder(pr, c("gene", "transcript", "aa_pos", "spp_1", "spp_2", "class", "type"))
        pr[]
    }
