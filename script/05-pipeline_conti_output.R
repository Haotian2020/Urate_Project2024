source("fn-ld_clump_local.R")
source("fn-ukbgwas-conti-output-format.R")

# for ukb full sample gwas output format ---------------------------------------
# urate gwas format and instruments --------------------------------------------

urateclean = ContiOutputFormatFunction(dat_name = "df_urate_clean_imputed.txt.gz",outcome_name = "urate_clean",
                                       mean = 308.834600400485, sd = 79.895501563483, samplesize=440466)
urate_tophits = ld_clump_local(urateclean)

# write files ------------------------------------------------------------------

write_tsv(urate_tophits, path = paste0(rdsf_personal, 'data/format_data/urate_clean_tophits.tsv'))
write.table(urateclean, file = paste0(rdsf_personal,"data/format_data/urate_clean_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

# sbp gwas format and instruments ----------------------------------------------

sbpclean = ContiOutputFormatFunction(dat_name = "df_sbp_clean_imputed.txt.gz",outcome_name = "sbp_clean",
                                     mean = 139.78564403065, sd = 19.5241937805439, samplesize=432099)
sbp_tophits = ld_clump_local(sbpclean)

# write files ------------------------------------------------------------------

write_tsv(sbp_tophits, path = paste0(rdsf_personal, 'data/format_data/sbp_clean_tophits.tsv'))
write.table(sbpclean, file = paste0(rdsf_personal,"data/format_data/sbp_clean_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

# dbp gwas format and instruments ----------------------------------------------

dbpclean = ContiOutputFormatFunction(dat_name = "df_dbp_clean_imputed.txt.gz",outcome_name = "dbp_clean",
                                     mean = 82.1392934230647, sd = 10.6322476240864, samplesize=432253 )
dbp_tophits = ld_clump_local(dbpclean)

# write files ------------------------------------------------------------------
write_tsv(dbp_tophits, path = paste0(rdsf_personal, 'data/format_data/dbp_clean_tophits.tsv'))
write.table(dbpclean, file = paste0(rdsf_personal,"data/format_data/dbp_clean_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)



# for ukb randomized sample gwas output format ---------------------------------
# splitted samples: urate gwas format and instruments --------------------------

urates1clean = ContiOutputFormatFunction(dat_name = "df_urate_clean_s1_imputed.txt.gz",outcome_name = "urate_clean_s1",
                                         mean = 308.886827182596, sd = 80.0756188911649, samplesize=220082)
urates1_tophits = ld_clump_local(urates1clean)

urates2clean = ContiOutputFormatFunction(dat_name = "df_urate_clean_s2_imputed.txt.gz",outcome_name = "urate_clean_s2",
                                         mean = 308.744211449421, sd = 79.6988890362597, samplesize=220081 )
urates2_tophits = ld_clump_local(urates2clean)

# write files ------------------------------------------------------------------

write_tsv(urates1_tophits, path = paste0(rdsf_personal, 'data/format_data/urate_clean_s1_tophits.tsv'))
write.table(urates1clean, file = paste0(rdsf_personal,"data/format_data/urate_clean_s1_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

write_tsv(urates2_tophits, path = paste0(rdsf_personal, 'data/format_data/urate_clean_s2_tophits.tsv'))
write.table(urates2clean, file = paste0(rdsf_personal,"data/format_data/urate_clean_s2_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

# splitted samples: sbp gwas format and instruments ----------------------------

sbps1clean = ContiOutputFormatFunction(dat_name = "df_sbp_clean_s1_imputed.txt.gz",outcome_name = "sbp_clean_s1",
                                       mean = 139.766111866107, sd = 19.4926971630032, samplesize=205299 )
sbps1_tophits = ld_clump_local(sbps1clean)

sbps2clean = ContiOutputFormatFunction(dat_name = "df_sbp_clean_s2_imputed.txt.gz",outcome_name = "sbp_clean_s2",
                                       mean = 139.777736161105, sd = 19.5385625835847, samplesize=205580 )
sbps2_tophits = ld_clump_local(sbps2clean)

# write files ------------------------------------------------------------------

write_tsv(sbps1_tophits, path = paste0(rdsf_personal, 'data/format_data/sbp_clean_s1_tophits.tsv'))
write.table(sbps1clean, file = paste0(rdsf_personal,"data/format_data/sbp_clean_s1_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

write_tsv(sbps2_tophits, path = paste0(rdsf_personal, 'data/format_data/sbp_clean_s2_tophits.tsv'))
write.table(sbps2clean, file = paste0(rdsf_personal,"data/format_data/sbp_clean_s2_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

# splitted samples: dbp gwas format and instruments ----------------------------

dbps1clean = ContiOutputFormatFunction(dat_name = "df_dbp_clean_s1_imputed.txt.gz",outcome_name = "dbp_clean_s1",
                                       mean = 82.1248745359047, sd = 10.6264048567287, samplesize=205238 )
dbps1_tophits = ld_clump_local(dbps1clean)

dbps2clean = ContiOutputFormatFunction(dat_name = "df_dbp_clean_s2_imputed.txt.gz",outcome_name = "dbp_clean_s2",
                                       mean = 82.1314181167036, sd = 10.5985958366052, samplesize=205512 )
dbps2_tophits = ld_clump_local(dbps2clean)

# write files ------------------------------------------------------------------

write_tsv(dbps1_tophits, path = paste0(rdsf_personal, 'data/format_data/dbp_clean_s1_tophits.tsv'))
write.table(dbps1clean, file = paste0(rdsf_personal,"data/format_data/dbp_clean_s1_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

write_tsv(dbps2_tophits, path = paste0(rdsf_personal, 'data/format_data/dbp_clean_s2_tophits.tsv'))
write.table(dbps2clean, file = paste0(rdsf_personal,"data/format_data/dbp_clean_s2_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)
