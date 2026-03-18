
library(data.table)

out_dir <- "/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/classification/"
out_csv <- paste0(out_dir, "shared_classA.csv")
out_rds <- paste0(out_dir, "shared_classA.rds")
if(!file.exists(out_csv)) file.create(out_csv)
if(!file.exists(out_rds)) file.create(out_rds)

rds_file <- paste0(out_dir, "shared_dt_test.rds")

shared_dt <- readRDS(rds_file)
message(nrow(shared_dt), " total variants in shared table")

### shared table rds headings
# chr	pos	ref	alt	
# variant.id_mel	af_mel	effect_order_mel	effect_mel	impact_mel	gene_mel	gene_id_mel	feature_type_mel	transcript_id_mel	biotype_mel	in_exon_mel	nt_change_mel	aa_change_mel	aa_pos_mel	aa_sub_mel	
# variant.id_sim	af_sim	effect_order_sim	effect_sim	impact_sim	gene_sim	gene_id_sim	feature_type_sim	transcript_id_sim	biotype_sim	in_exon_sim	nt_change_sim	aa_change_sim	aa_pos_sim	aa_sub_sim

### get columns
af_mel <- "af_mel"
af_sim <- "af_sim"
aa_mel <- "aa_sub_mel"
aa_sim <- "aa_sub_sim"

### MAKE A COPY of shared dt where 
# variant exists for both species
#  and amino acid is the same 
shared_classA <- shared_dt[
    !is.na(get(af_mel)) & !is.na(get(af_sim)) & get(aa_mel) == get(aa_sim)
]
message(nrow(shared_classA), " shared variants with same pos, nt, and aa")

saveRDS(shared_classA, out_rds)
message("RDS saved to ", out_rds)

if(nrow(shared_classA) < 600){
    fwrite(shared_classA, out_csv)
    message("CSV saved to ", out_csv)
}
