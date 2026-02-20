library(SeqArray)
library(data.table)

shared_table <- readRDS("/scratch/ejy4bu/drosophila/gds_analysis/snp_datatables/shared_snp_dt_table.rds")
head(shared_table)
dim(shared_table)
str(shared_table)

# any NAs in key columns?
shared_table[, lapply(.SD, function(x) sum(is.na(x)))]

# what effects are present?
shared_table[, .N, by = effect][order(-N)]



# > head(shared_table)
# table)   variant.id    chr   pos    ref    alt     af_mel    af_sim
#         <int> <char> <int> <char> <char>      <num>     <num>
# 1:       1162     2L 10247      C      G 0.00236771 0.0045045
# 2:       1162     2L 10247      C      G 0.00236771 0.0045045
# 3:       1162     2L 10247      C      G 0.00236771 0.0045045
# 4:       1162     2L 10247      C      G 0.00236771 0.0045045
# 5:       1162     2L 10247      C      G 0.00236771 0.0045045
# 6:       1162     2L 10247      C      G 0.00236771 0.0045045
#                 effect   impact   gene     gene_id feature_type transcript_id
#                 <char>   <char> <char>      <char>       <char>        <char>
# 1: 3_prime_UTR_variant MODIFIER l(2)gl FBgn0002121   transcript   FBtr0078170
# 2: 3_prime_UTR_variant MODIFIER l(2)gl FBgn0002121   transcript   FBtr0078171
# 3: 3_prime_UTR_variant MODIFIER l(2)gl FBgn0002121   transcript   FBtr0078166
# 4: 3_prime_UTR_variant MODIFIER l(2)gl FBgn0002121   transcript   FBtr0078167
# 5: 3_prime_UTR_variant MODIFIER l(2)gl FBgn0002121   transcript   FBtr0078168
# 6: 3_prime_UTR_variant MODIFIER l(2)gl FBgn0002121   transcript   FBtr0078169
#           biotype in_exon nt_change aa_change aa_pos
#            <char>  <char>    <char>    <char> <char>
# 1: protein_coding     9/9 c.*968G>C
# 2: protein_coding   10/10 c.*968G>C
# 3: protein_coding     9/9 c.*968G>C
# 4: protein_coding     9/9 c.*968G>C
# 5: protein_coding     9/9 c.*968G>C
# 6: protein_coding   10/10 c.*968G>C



# > dim(shared_table)
# [1] 7603036      18



# > str(shared_table)
# Classes ‘data.table’ and 'data.frame':  7603036 obs. of  18 variables:
#  $ variant.id   : int  1162 1162 1162 1162 1162 1162 1162 1162 1162 1162 ...
#  $ chr          : chr  "2L" "2L" "2L" "2L" ...
#  $ pos          : int  10247 10247 10247 10247 10247 10247 10247 10247 10247 10247 ...
#  $ ref          : chr  "C" "C" "C" "C" ...
#  $ alt          : chr  "G" "G" "G" "G" ...
#  $ af_mel       : num  0.00237 0.00237 0.00237 0.00237 0.00237 ...
#  $ af_sim       : num  0.0045 0.0045 0.0045 0.0045 0.0045 ...
#  $ effect       : chr  "3_prime_UTR_variant" "3_prime_UTR_variant" "3_prime_UTR_variant" "3_prime_UTR_variant" ...
#  $ impact       : chr  "MODIFIER" "MODIFIER" "MODIFIER" "MODIFIER" ...
#  $ gene         : chr  "l(2)gl" "l(2)gl" "l(2)gl" "l(2)gl" ...
#  $ gene_id      : chr  "FBgn0002121" "FBgn0002121" "FBgn0002121" "FBgn0002121" ...
#  $ feature_type : chr  "transcript" "transcript" "transcript" "transcript" ...
#  $ transcript_id: chr  "FBtr0078170" "FBtr0078171" "FBtr0078166" "FBtr0078167" ...
#  $ biotype      : chr  "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
#  $ in_exon      : chr  "9/9" "10/10" "9/9" "9/9" ...
#  $ nt_change    : chr  "c.*968G>C" "c.*968G>C" "c.*968G>C" "c.*968G>C" ...
#  $ aa_change    : chr  "" "" "" "" ...
#  $ aa_pos       : chr  "" "" "" "" ...
#  - attr(*, ".internal.selfref")=<externalptr>
# >


# > # any NAs in key columns?
# > shared_table[, lapply(.SD, function(x) sum(is.na(x)))]

#    variant.id   chr   pos   ref   alt af_mel af_sim effect impact  gene gene_id
#         <int> <int> <int> <int> <int>  <int>  <int>  <int>  <int> <int>   <int>
# 1:          0     0     0     0     0      0      0      0      0     0       0
#    feature_type transcript_id biotype in_exon nt_change aa_change aa_pos
#           <int>         <int>   <int>   <int>     <int>     <int>  <int>
# 1:            0             0       0       0         0         0      0
# >


# > # what effects are present?
# > shared_table[, .N, by = effect][order(-N)]
#                                                       effect       N
#                                                       <char>   <int>
#  1:                                           intron_variant 2304292
#  2:                                    upstream_gene_variant 1993091
#  3:                                  downstream_gene_variant 1982432
#  4:                                       synonymous_variant  477113
#  5:                                        intergenic_region  333147
#  6:                                         missense_variant  223109
#  7:                                      3_prime_UTR_variant  117090
#  8:                                      5_prime_UTR_variant   77493
#  9:                     splice_region_variant&intron_variant   33034
# 10:                       non_coding_transcript_exon_variant   30787
# 11:                                              stop_gained   12744
# 12:           5_prime_UTR_premature_start_codon_gain_variant    9438
# 13:                 splice_region_variant&synonymous_variant    3315
# 14:                   missense_variant&splice_region_variant    2030
# 15:                                    splice_region_variant    1570
# 16:                      splice_donor_variant&intron_variant     651
# 17:                   splice_acceptor_variant&intron_variant     374
# 18:                                                stop_lost     362
# 19: splice_region_variant&non_coding_transcript_exon_variant     300
# 20:                                    stop_retained_variant     212
# 21:                        stop_gained&splice_region_variant     212
# 22:                                               start_lost     137
# 23:                                  initiator_codon_variant      84
# 24:              splice_region_variant&stop_retained_variant       9
# 25:                          stop_lost&splice_region_variant       5
# 26:            initiator_codon_variant&splice_region_variant       4
# 27:                         start_lost&splice_region_variant       1
#                                                       effect       N
# >
# >
