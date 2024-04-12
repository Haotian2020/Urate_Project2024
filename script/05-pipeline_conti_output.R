# Script for formatting GWAS summary level from UKB pipeline output of continuous traits

source("fn-ld_clump_local.R")
source("fn-ukbgwas-conti-output-format.R")

# for ukb full sample gwas output format ---------------------------------------
# urate gwas format and instruments --------------------------------------------

urateclean = ContiOutputFormatFunction(dat_name = "df_urate_clean_imputed.txt.gz",outcome_name = "Urate (UKB)",
                                       mean = 308.834600400485, sd = 79.895501563483, samplesize=440466)
urate_tophits = ld_clump_local(urateclean)

# write files ------------------------------------------------------------------

write_tsv(urate_tophits, file = paste0(rdsf_personal, 'data/format_data/urate_clean_tophits.tsv'))
data.table::fwrite(urateclean,paste0(rdsf_personal,"data/format_data/urate_clean_GWAS_tidy_outcome.csv"))

# sbp gwas format and instruments ----------------------------------------------

sbpclean = ContiOutputFormatFunction(dat_name = "df_sbp_clean_imputed.txt.gz",outcome_name = "SBP (UKB)",
                                     mean = 139.78564403065, sd = 19.5241937805439, samplesize=432099)
sbp_tophits = ld_clump_local(sbpclean)

# write files ------------------------------------------------------------------

write_tsv(sbp_tophits, file = paste0(rdsf_personal, 'data/format_data/sbp_clean_tophits.tsv'))
data.table::fwrite(sbpclean,paste0(rdsf_personal,"data/format_data/sbp_clean_GWAS_tidy_outcome.csv"))

# dbp gwas format and instruments ----------------------------------------------

dbpclean = ContiOutputFormatFunction(dat_name = "df_dbp_clean_imputed.txt.gz",outcome_name = "DBP (UKB)",
                                     mean = 82.1392934230647, sd = 10.6322476240864, samplesize=432253 )
dbp_tophits = ld_clump_local(dbpclean)

# write files ------------------------------------------------------------------
write_tsv(dbp_tophits, file = paste0(rdsf_personal, 'data/format_data/dbp_clean_tophits.tsv'))
data.table::fwrite(dbpclean,paste0(rdsf_personal,"data/format_data/dbp_clean_GWAS_tidy_outcome.csv"))

# for ukb randomized sample gwas output format ---------------------------------
# splitted samples: urate gwas format and instruments --------------------------

urates1clean = ContiOutputFormatFunction(dat_name = "df_urate_clean_s1_imputed.txt.gz",outcome_name = "Urate (UKB s1)",
                                         mean = 308.886827182596, sd = 80.0756188911649, samplesize=220082)
urates1_tophits = ld_clump_local(urates1clean)

urates2clean = ContiOutputFormatFunction(dat_name = "df_urate_clean_s2_imputed.txt.gz",outcome_name = "Urate (UKB s2)",
                                         mean = 308.744211449421, sd = 79.6988890362597, samplesize=220081 )
urates2_tophits = ld_clump_local(urates2clean)

# write files ------------------------------------------------------------------

write_tsv(urates1_tophits, file = paste0(rdsf_personal, 'data/format_data/urate_clean_s1_tophits.tsv'))
data.table::fwrite(urates1clean,paste0(rdsf_personal,"data/format_data/urate_clean_s1_GWAS_tidy_outcome.csv"))

write_tsv(urates2_tophits, file = paste0(rdsf_personal, 'data/format_data/urate_clean_s2_tophits.tsv'))
data.table::fwrite(urates2clean,paste0(rdsf_personal,"data/format_data/urate_clean_s2_GWAS_tidy_outcome.csv"))

# splitted samples: sbp gwas format and instruments ----------------------------

sbps1clean = ContiOutputFormatFunction(dat_name = "df_sbp_clean_s1_imputed.txt.gz",outcome_name = "SBP (UKB s1)",
                                       mean = 139.766111866107, sd = 19.4926971630032, samplesize=205299 )
sbps1_tophits = ld_clump_local(sbps1clean)

sbps2clean = ContiOutputFormatFunction(dat_name = "df_sbp_clean_s2_imputed.txt.gz",outcome_name = "SBP (UKB s2)",
                                       mean = 139.777736161105, sd = 19.5385625835847, samplesize=205580 )
sbps2_tophits = ld_clump_local(sbps2clean)

# write files ------------------------------------------------------------------

write_tsv(sbps1_tophits, file = paste0(rdsf_personal, 'data/format_data/sbp_clean_s1_tophits.tsv'))
data.table::fwrite(sbps1clean,paste0(rdsf_personal,"data/format_data/sbp_clean_s1_GWAS_tidy_outcome.csv"))

write_tsv(sbps2_tophits, file = paste0(rdsf_personal, 'data/format_data/sbp_clean_s2_tophits.tsv'))
data.table::fwrite(sbps2clean,paste0(rdsf_personal,"data/format_data/sbp_clean_s2_GWAS_tidy_outcome.csv"))

# splitted samples: dbp gwas format and instruments ----------------------------

dbps1clean = ContiOutputFormatFunction(dat_name = "df_dbp_clean_s1_imputed.txt.gz",outcome_name = "DBP (UKB s1)",
                                       mean = 82.1248745359047, sd = 10.6264048567287, samplesize=205238 )
dbps1_tophits = ld_clump_local(dbps1clean)

dbps2clean = ContiOutputFormatFunction(dat_name = "df_dbp_clean_s2_imputed.txt.gz",outcome_name = "DBP (UKB s2)",
                                       mean = 82.1314181167036, sd = 10.5985958366052, samplesize=205512 )
dbps2_tophits = ld_clump_local(dbps2clean)

# write files ------------------------------------------------------------------

write_tsv(dbps1_tophits, file = paste0(rdsf_personal, 'data/format_data/dbp_clean_s1_tophits.tsv'))
data.table::fwrite(dbps1clean,paste0(rdsf_personal,"data/format_data/dbp_clean_s1_GWAS_tidy_outcome.csv"))

write_tsv(dbps2_tophits, file = paste0(rdsf_personal, 'data/format_data/dbp_clean_s2_tophits.tsv'))
data.table::fwrite(dbps2clean,paste0(rdsf_personal,"data/format_data/dbp_clean_s2_GWAS_tidy_outcome.csv"))
