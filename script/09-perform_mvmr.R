source("fn-mvmr.R")
source("fn-mvmr_meta.R")

# external urate and bp (UKB) on egfr ckdgen -----------------------------------

mvmr_exurate_sbp_egfr = MVMR_function("exurate_sd","sbp_clean","egfr_sd")

print(mvmr_instruments_strength("exurate_sd","sbp_clean","egfr_sd"))

mvmr_exurate_dbp_egfr = MVMR_function("exurate_sd","dbp_clean","egfr_sd")

print(mvmr_instruments_strength("exurate_sd","dbp_clean","egfr_sd"))

# format names -----------------------------------------------------------------

mvmr_exurate_sbp_egfr$result$method = 'MVMR'
exurate_sbp_egfr_mvmr = data.frame(mvmr_exurate_sbp_egfr$result)
exurate_sbp_egfr_mvmr$outcome = "eGFR (CKDGen)"

# format names -----------------------------------------------------------------
mvmr_exurate_dbp_egfr$result$method = 'MVMR'
exurate_dbp_egfr_mvmr = data.frame(mvmr_exurate_dbp_egfr$result)
exurate_dbp_egfr_mvmr$outcome = "eGFR (CKDGen)"

# write results ----------------------------------------------------------------
# mvmr results need to be save seperately to distinguish the duplicated exposure effects

write.table(exurate_sbp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_sbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(exurate_dbp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_dbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)

# urate and bp from UKB on egfr ckdgen -----------------------------------------

mvmr_urates1_sbps2_egfr = MVMR_function("urate_clean_s1","sbp_clean_s2","egfr_sd")
mvmr_urates1_sbps2_egfr$method = 'MVMR'
mvmr_urates2_sbps1_egfr = MVMR_function("urate_clean_s2","sbp_clean_s1","egfr_sd")
mvmr_urates2_sbps1_egfr$method = 'MVMR'

uratesbp_egfr_mvmr = format_mvmrmeta(mvmr_urates1_sbps2_egfr,mvmr_urates2_sbps1_egfr,"Urate","SBP")

print(mvmr_instruments_strength("urate_clean_s1","sbp_clean_s2","egfr_sd"))

print(mvmr_instruments_strength("urate_clean_s2","sbp_clean_s1","egfr_sd"))

# format names -----------------------------------------------------------------

uratesbp_egfr_mvmr["exposure"][uratesbp_egfr_mvmr["exposure"] == "Urate (UKB s1) + Urate (UKB s2)"] <- "Urate (UKB Meta)"
uratesbp_egfr_mvmr["exposure"][uratesbp_egfr_mvmr["exposure"] == "SBP (UKB s2) + SBP (UKB s1)"] <- "SBP (UKB Meta)"
uratesbp_egfr_mvmr$outcome = "eGFR (CKDGen)"

# 
mvmr_urates1_dbps2_egfr = MVMR_function("urate_clean_s1","dbp_clean_s2","egfr_sd")
mvmr_urates1_dbps2_egfr$method = 'MVMR'
mvmr_urates2_dbps1_egfr = MVMR_function("urate_clean_s2","dbp_clean_s1","egfr_sd")
mvmr_urates2_dbps1_egfr$method = 'MVMR'

uratedbp_egfr_mvmr = format_mvmrmeta(mvmr_urates1_dbps2_egfr,mvmr_urates2_dbps1_egfr,"Urate","DBP")

mvmr_instruments_strength("urate_clean_s1","dbp_clean_s2","egfr_sd")

mvmr_instruments_strength("urate_clean_s2","dbp_clean_s1","egfr_sd")

# format names -----------------------------------------------------------------

uratedbp_egfr_mvmr["exposure"][uratedbp_egfr_mvmr["exposure"] == "Urate (UKB s1) + Urate (UKB s2)"] <- "Urate (UKB Meta)"
uratedbp_egfr_mvmr["exposure"][uratedbp_egfr_mvmr["exposure"] == "DBP (UKB s2) + DBP (UKB s1)"] <- "DBP (UKB Meta)"
uratedbp_egfr_mvmr$outcome = "eGFR (CKDGen)"

# write results ----------------------------------------------------------------

write.table(uratesbp_egfr_mvmr, file = paste0(rdsf_personal,"results/uratesbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(uratedbp_egfr_mvmr, file = paste0(rdsf_personal,"results/uratedbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)

# addtional test with bmi ------------------------------------------------------

mvmr_sbp_bmi_egfr = MVMR_function("sbp_clean","ieu-a-2","egfr_sd")

mvmr_dbp_bmi_egfr = MVMR_function("dbp_clean","ieu-a-2","egfr_sd")
