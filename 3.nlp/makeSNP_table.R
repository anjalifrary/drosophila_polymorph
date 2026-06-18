library(SeqArray)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(16)


genofile <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.sim.all.SNAPE.001.50.20Nov2025_sim.norep.ann.dmel6.eff.gds")
# genofile <- seqOpen("/scratch/ejy4bu/drosophila/gds_files/dest.PoolSeq.SNAPE.001.50.03Dec2024_DACtest.norep.ann.eff.gds")

codon_mutation_matrix <- function() {
# Standard genetic code
genetic_code <- c(
    TTT="F", TTC="F", TTA="L", TTG="L",
    TCT="S", TCC="S", TCA="S", TCG="S",
    TAT="Y", TAC="Y", TAA="*", TAG="*",
    TGT="C", TGC="C", TGA="*", TGG="W",
    CTT="L", CTC="L", CTA="L", CTG="L",
    CCT="P", CCC="P", CCA="P", CCG="P",
    CAT="H", CAC="H", CAA="Q", CAG="Q",
    CGT="R", CGC="R", CGA="R", CGG="R",
    ATT="I", ATC="I", ATA="I", ATG="M",
    ACT="T", ACC="T", ACA="T", ACG="T",
    AAT="N", AAC="N", AAA="K", AAG="K",
    AGT="S", AGC="S", AGA="R", AGG="R",
    GTT="V", GTC="V", GTA="V", GTG="V",
    GCT="A", GCC="A", GCA="A", GCG="A",
    GAT="D", GAC="D", GAA="E", GAG="E",
    GGT="G", GGC="G", GGA="G", GGG="G"
)

codons <- names(genetic_code)
bases <- c("A","C","G","T")

# --- Generate all single mutants (vectorized) ---
# Create 3 mutation positions × 3 alternative bases = 9 mutations per codon
pos <- rep(1:3, each = 3)
alt_base_idx <- rep(1:4, times = 3)

# For each codon, generate all substitutions
mutate_codon <- function(codon) {
    chars <- strsplit(codon, "")[[1]]

    # For each position, generate 3 alternative bases
    muts <- unlist(lapply(1:3, function(i) {
    setdiff(bases, chars[i]) |>
        sapply(function(b) {
        tmp <- chars
        tmp[i] <- b
        paste(tmp, collapse = "")
        })
    }))

    return(muts)
}

# Apply across all codons (vectorized via sapply)
mutant_matrix <- sapply(codons, mutate_codon)  # 9 x 64 matrix

# Flatten
original_codons <- rep(codons, each = 9)
mutant_codons <- as.vector(mutant_matrix)

# Map to amino acids
aa_original <- genetic_code[original_codons]
aa_mutant <- genetic_code[mutant_codons]

# Classify mutations
synonymous <- aa_original == aa_mutant
nonsense <- aa_mutant == "*"
nonsynonymous <- !(synonymous | nonsense)

# Build data frame
df <- data.frame(
    codon = original_codons,
    synonymous = synonymous,
    nonsynonymous = nonsynonymous,
    nonsense = nonsense
)

# Aggregate per codon
per_codon <- aggregate(. ~ codon, data = df, mean)

# Genome-wide expectation (sense codons only)
sense <- per_codon[genetic_code[per_codon$codon] != "*", ]

genome_avg <- colMeans(sense[, c("synonymous","nonsynonymous","nonsense")])

return(list(
    per_codon = per_codon,
    genome_average = genome_avg
))
}


res <- codon_mutation_matrix()

### common SNP.dt
message("Making basic SNP table")
seqResetFilter(genofile)
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                    pos=seqGetData(genofile, "position"),
                    nAlleles=seqGetData(genofile, "$num_allele"),
                    id=seqGetData(genofile, "variant.id"))

snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
seqSetFilter(genofile, snp.dt$id)
snp.dt[,bin:=floor(id/10000)]
table(snp.dt$chr)

### get annotations
message("Getting annotations")
ann <- foreach(i=unique(snp.dt$bin), .combine="rbind")%dopar%{
# i <- 100
message(paste(i, max(snp.dt$bin), sep=" / "))
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id=snp.dt[bin==i]$id)

tmp <- seqGetData(genofile, "annotation/info/EFF")
len1 <- tmp$length
len2 <- tmp$data

snp.dt1 <- data.table(len=rep(len1, times=len1),
                        ann=len2,
                        id=rep(snp.dt[bin==i]$id, times=len1))

snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\(")[[1]]]
snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[6]]]
snp.dt1[,codons:=tstrsplit(snp.dt1$ann,"\\|")[[3]]]
snp.dt1[,aas:=tstrsplit(tstrsplit(snp.dt1$ann,"\\|")[[4]], "/")[[1]]]
snp.dt1[,transcript:=tstrsplit(snp.dt1$ann,"\\|")[[9]]]

snp.dt1 <- snp.dt1[class%in%c("missense_variant", "synonymous_variant")]


snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=","), codons=paste(codons, collapse=","), aas=paste(aas, collapse=","), transcript=paste(transcript, collapse=",")),
                        list(variant.id=id)]

snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]
snp.dt1.an[,codons:=tstrsplit(snp.dt1.an$codons,"\\,")[[1]]]
snp.dt1.an[,transcript:=tstrsplit(snp.dt1.an$transcript,"\\,")[[1]]]

snp.dt1.an[,ref_codon:=tstrsplit(codons, "/")[[1]]]
snp.dt1.an[,alt_codon:=tstrsplit(codons, "/")[[2]]]
snp.dt1.an[,aas:=tstrsplit(snp.dt1.an$aas,"\\,")[[1]]]
snp.dt1.an[,ref_aa:=sub("^p\\.([A-Za-z]{3})\\d+[A-Za-z]{3}$", "\\1", snp.dt1.an$aas)]
snp.dt1.an[,alt_aa:=sub("^p\\.[A-Za-z]{3}\\d+([A-Za-z]{3})$", "\\1", snp.dt1.an$aas)]
snp.dt1.an[,aa_pos:=as.numeric(sub("^p\\.[A-Za-z]{3}(\\d+)[A-Za-z]{3}$", "\\1", snp.dt1.an$aas))]
snp.dt1.an[,codon:=toupper(ref_codon)]


eff <- merge(snp.dt1.an, res$per_codon, by="codon")
eff[,c("variant.id", "col", "gene", "transcript", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "aa_pos", "synonymous", "nonsynonymous", "nonsense")]

return(eff)
}


snp.dt <- merge(snp.dt, ann, by.x="id", by.y="variant.id")



### get some basic filtering info
message("Basic filtering")
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id=snp.dt$id)

try(snp.dt[,filter:=seqGetData(genofile, "annotation/filter")])

snp.dt[,ref:=seqGetData(genofile, "$ref")]
snp.dt[,alt:=seqGetData(genofile, "$alt")]

setkey(snp.dt, "ref")
snp.dt[,b:=0]
snp.dt[J(c("A", "C", "T", "G")), b:=1]
setkey(snp.dt, "alt")
snp.dt[J(c("A", "C", "T", "G")), b:=b+1]
snp.dt[b==2]
snp.dt[,use:=F]
snp.dt[b==2, use:=T]
table(snp.dt$use)
table(snp.dt$filter)

### try to get the standardized quality info
QUAL_columns <- c("ExcessHet","FS","InbreedingCoeff", "ICB", "HOB")
#QUAL_columns <- c("foo", "bar")

qual.tmp <- foreach(qq=QUAL_columns, .errorhandling="remove", .combine="cbind")%do%{
## qq <- QUAL_columns[1]
tmp <- data.table(qq=seqGetData(genofile, paste("annotation/info", qq, sep="/")))
setnames(tmp, "qq", qq)
return(tmp)
}
#str(qual.tmp)
#str(snp.dt)
try(snp.dt <- cbind(snp.dt, qual.tmp))
#str(snp.dt)

### flag adjacent SNPs in the same codon
snp.dt.codon <- snp.dt[,list(mutations_per_codon=.N), list(aa_pos, transcript)]
setkey(snp.dt.codon, aa_pos, transcript)
setkey(snp.dt, aa_pos, transcript)

snp.dt <- merge(snp.dt.codon, snp.dt)
setkey(snp.dt, chr, pos)
snp.dt[mutations_per_codon==3]
prop.table(table(snp.dt.codon$N))
