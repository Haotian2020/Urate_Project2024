# Script for formatting GWAS summary level from UKB pipeline output of continuous traits

source("fn-ld_clump_local.R")
source("fn-ukbgwas-conti-output-format.R")

# for ukb full sample gwas output format ---------------------------------------
# urate gwas format and instruments --------------------------------------------

urateclean = ContiOutputFormatFunction(dat_name = "urate_imputed.txt.gz",outcome_name = "Urate (UKB)",
                                       mean = 308.834600400485, sd = 79.895501563483, samplesize=440466)
urate_tophits = ld_clump_local(urateclean)

# write files ------------------------------------------------------------------

write_tsv(urate_tophits, file = paste0(rdsf_personal, 'data/format_data/urate_tophits.tsv'))
data.table::fwrite(urateclean,paste0(rdsf_personal,"data/format_data/urate_GWAS_tidy_outcome.csv"))

# sbp gwas format and instruments ----------------------------------------------

sbpclean = ContiOutputFormatFunction(dat_name = "sbp_avg_imputed.txt.gz",outcome_name = "SBP (UKB)",
                                     mean = 137.861762634246, sd = 18.4833994911214, samplesize = 435859)
sbp_tophits = ld_clump_local(sbpclean)

# write files ------------------------------------------------------------------

write_tsv(sbp_tophits, file = paste0(rdsf_personal, 'data/format_data/sbp_tophits.tsv'))
data.table::fwrite(sbpclean,paste0(rdsf_personal,"data/format_data/sbp_GWAS_tidy_outcome.csv"))

# dbp gwas format and instruments ----------------------------------------------

dbpclean = ContiOutputFormatFunction(dat_name = "dbp_avg_imputed.txt.gz",outcome_name = "DBP (UKB)",
                                     mean = 82.1555713017935, sd = 10.0886399141079, samplesize = 436083)
dbp_tophits = ld_clump_local(dbpclean)

# write files ------------------------------------------------------------------
write_tsv(dbp_tophits, file = paste0(rdsf_personal, 'data/format_data/dbp_tophits.tsv'))
data.table::fwrite(dbpclean,paste0(rdsf_personal,"data/format_data/dbp_GWAS_tidy_outcome.csv"))

# pp gwas format and instruments -----------------------------------------------

ppclean = ContiOutputFormatFunction(dat_name = "pp_avg_imputed.txt.gz",outcome_name = "PP (UKB)",
                                     mean = 55.6526896881128, sd = 13.4170610930985, samplesize = 435478 )
pp_tophits = ld_clump_local(ppclean)

# write files ------------------------------------------------------------------
write_tsv(pp_tophits, file = paste0(rdsf_personal, 'data/format_data/pp_tophits.tsv'))
data.table::fwrite(ppclean, paste0(rdsf_personal,"data/format_data/pp_GWAS_tidy_outcome.csv"))

# for ukb randomized sample gwas output format ---------------------------------
# splitted samples: urate gwas format and instruments --------------------------

urates1clean = ContiOutputFormatFunction(dat_name = "urate_s1_imputed.txt.gz",outcome_name = "Urate (UKB s1)",
                                         mean = 308.7965, sd = 79.78696, samplesize = 206924)
urates1_tophits = ld_clump_local(urates1clean)

urates2clean = ContiOutputFormatFunction(dat_name = "urate_s2_imputed.txt.gz",outcome_name = "Urate (UKB s2)",
                                         mean = 308.9642, sd = 79.90702, samplesize = 206923 )
urates2_tophits = ld_clump_local(urates2clean)

# write files ------------------------------------------------------------------

write_tsv(urates1_tophits, file = paste0(rdsf_personal, 'data/format_data/urate_s1_tophits.tsv'))
data.table::fwrite(urates1clean,paste0(rdsf_personal,"data/format_data/urate_s1_GWAS_tidy_outcome.csv"))

write_tsv(urates2_tophits, file = paste0(rdsf_personal, 'data/format_data/urate_s2_tophits.tsv'))
data.table::fwrite(urates2clean,paste0(rdsf_personal,"data/format_data/urate_s2_GWAS_tidy_outcome.csv"))

# splitted samples: sbp gwas format and instruments ----------------------------

sbps1clean = ContiOutputFormatFunction(dat_name = "sbp_avg_s1_imputed.txt.gz",outcome_name = "SBP (UKB s1)",
                                       mean = 137.744, sd = 18.31552, samplesize = 206924 )
sbps1_tophits = ld_clump_local(sbps1clean)

sbps2clean = ContiOutputFormatFunction(dat_name = "sbp_avg_s2_imputed.txt.gz",outcome_name = "SBP (UKB s2)",
                                       mean = 137.7833, sd = 18.3726, samplesize = 206923 )
sbps2_tophits = ld_clump_local(sbps2clean)

# write files ------------------------------------------------------------------

write_tsv(sbps1_tophits, file = paste0(rdsf_personal, 'data/format_data/sbp_s1_tophits.tsv'))
data.table::fwrite(sbps1clean,paste0(rdsf_personal,"data/format_data/sbp_s1_GWAS_tidy_outcome.csv"))

write_tsv(sbps2_tophits, file = paste0(rdsf_personal, 'data/format_data/sbp_s2_tophits.tsv'))
data.table::fwrite(sbps2clean,paste0(rdsf_personal,"data/format_data/sbp_s2_GWAS_tidy_outcome.csv"))

# splitted samples: dbp gwas format and instruments ----------------------------

dbps1clean = ContiOutputFormatFunction(dat_name = "dbp_avg_s1_imputed.txt.gz",outcome_name = "DBP (UKB s1)",
                                       mean = 82.12938, sd = 10.04691, samplesize = 206924)
dbps1_tophits = ld_clump_local(dbps1clean)

dbps2clean = ContiOutputFormatFunction(dat_name = "dbp_avg_s2_imputed.txt.gz",outcome_name = "DBP (UKB s2)",
                                       mean = 82.15788, sd = 10.08448, samplesize = 206923 )
dbps2_tophits = ld_clump_local(dbps2clean)

# write files ------------------------------------------------------------------

write_tsv(dbps1_tophits, file = paste0(rdsf_personal, 'data/format_data/dbp_s1_tophits.tsv'))
data.table::fwrite(dbps1clean,paste0(rdsf_personal,"data/format_data/dbp_s1_GWAS_tidy_outcome.csv"))

write_tsv(dbps2_tophits, file = paste0(rdsf_personal, 'data/format_data/dbp_s2_tophits.tsv'))
data.table::fwrite(dbps2clean,paste0(rdsf_personal,"data/format_data/dbp_s2_GWAS_tidy_outcome.csv"))


# splitted samples: pp gwas format and instruments ----------------------------

pps1clean = ContiOutputFormatFunction(dat_name = "pp_avg_s1_imputed.txt.gz",outcome_name = "PP (UKB s1)",
                                       mean = 55.6146, sd = 13.36189, samplesize = 206924 )
pps1_tophits = ld_clump_local(pps1clean)

pps2clean = ContiOutputFormatFunction(dat_name = "pp_avg_s2_imputed.txt.gz",outcome_name = "PP (UKB s2)",
                                       mean = 55.62538, sd = 13.40611, samplesize = 206923 )
pps2_tophits = ld_clump_local(pps2clean)

# write files ------------------------------------------------------------------

write_tsv(pps1_tophits, file = paste0(rdsf_personal, 'data/format_data/pp_s1_tophits.tsv'))
data.table::fwrite(pps1clean,paste0(rdsf_personal,"data/format_data/pp_s1_GWAS_tidy_outcome.csv"))

write_tsv(pps2_tophits, file = paste0(rdsf_personal, 'data/format_data/pp_s2_tophits.tsv'))
data.table::fwrite(pps2clean,paste0(rdsf_personal,"data/format_data/pp_s2_GWAS_tidy_outcome.csv"))
