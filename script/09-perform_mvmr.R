# Script for performing MVMR analyses

source("fn-mvmr.R")
source("fn-mvmr_meta.R")

# urate (CKDgen) and bp (UKB) on egfr (CKDGen2019) -----------------------------

mvmr_exurate_sbp_egfr = MVMR_function("exurate_sd","sbp","egfr_sd")

print(mvmr_exurate_sbp_egfr)
# urate sbp
# F-statistic  36.60349  44.97061

mvmr_exurate_dbp_egfr = MVMR_function("exurate_sd","dbp","egfr_sd")

print(mvmr_exurate_dbp_egfr)
# dbp urate
# F-statistic  42.15317  37.80472

mvmr_exurate_pp_egfr = MVMR_function("exurate_sd","pp","egfr_sd")

print(mvmr_exurate_pp_egfr)
# pp urate
# F-statistic  56.69814   32.6721

# urate and bmi on egfr (CKDGen2019) -------------------------------------------

mvmr_exurate_bmi_egfr = MVMR_function("exurate_sd","ukb-b-19953","egfr_sd")

mvmr_urate_bmi_egfr = MVMR_function("urate","ukb-b-19953","egfr_sd")

bmi_egfr = uvmr("ukb-b-19953", "egfr_sd")

bmi_egfr_cys = uvmr("ukb-b-19953", "egfr_cys_sd")
# 2           Weighted median  399 -0.3855276 0.02393745  2.329302e-58 Steiger
# 3 Inverse variance weighted  399 -0.3776102 0.01562088 4.235237e-129 Steiger
# 4               Simple mode  399 -0.3610444 0.07446347  1.786129e-06 Steiger
# 5             Weighted mode  399 -0.3934113 0.04105014  1.036006e-19 Steiger

# id.exposure id.outcome       outcome                                exposure
# 1 ukb-b-19953     xzp985 eGFR cys 2021 Body mass index (BMI) || id:ukb-b-19953
# egger_intercept           se      pval      Isq    type
# 1    0.0002807193 0.0007600161 0.7120566 0.984295 Steiger

# id.exposure id.outcome       outcome                                exposure
# 1 ukb-b-19953     xzp985 eGFR cys 2021 Body mass index (BMI) || id:ukb-b-19953
# 2 ukb-b-19953     xzp985 eGFR cys 2021 Body mass index (BMI) || id:ukb-b-19953
# method        Q Q_df       Q_pval    type
# 1                  MR Egger 738.2823  397 7.689325e-23 Steiger
# 2 Inverse variance weighted 738.5360  398 9.917593e-23 Steiger

MVMR_function("exurate_sd","ukb-b-19953","egfr_cys_sd")
# $result
# id.exposure                                exposure id.outcome     outcome
# 1      LrjzHo                          Urate (CKDGen)     gnf0Px egfr_cys_sd
# 2 ukb-b-19953 Body mass index (BMI) || id:ukb-b-19953     gnf0Px egfr_cys_sd
# nsnp          b         se         pval method
# 1   55 -0.1878733 0.03974361 2.277131e-06   MVMR
# 2  359 -0.3906078 0.03425483 4.039597e-30   MVMR

MVMR_function("urate","ukb-b-19953","egfr_cys_sd")
# $result
# id.exposure                                exposure id.outcome     outcome
# 1      jQPGtw                             Urate (UKB)     dWgHJ5 egfr_cys_sd
# 2 ukb-b-19953 Body mass index (BMI) || id:ukb-b-19953     dWgHJ5 egfr_cys_sd
# nsnp          b         se         pval method
# 1  194 -0.3529767 0.04342505 4.349519e-16   MVMR
# 2  270 -0.3319847 0.04601081 5.379038e-13   MVMR

# format names -----------------------------------------------------------------

exurate_sbp_egfr_mvmr = data.frame(mvmr_exurate_sbp_egfr$result)
exurate_sbp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

exurate_dbp_egfr_mvmr = data.frame(mvmr_exurate_dbp_egfr$result)
exurate_dbp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

exurate_pp_egfr_mvmr = data.frame(mvmr_exurate_pp_egfr$result)
exurate_pp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

exurate_bmi_egfr_mvmr = data.frame(mvmr_exurate_bmi_egfr$result)
exurate_bmi_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

# write results ----------------------------------------------------------------
# mvmr results need to be save separately to distinguish the duplicated exposure effects

write.table(exurate_sbp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_sbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(exurate_dbp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_dbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(exurate_pp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_pp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(exurate_bmi_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_bmi_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)

# urate and bp from UKB on egfr ckdgen -----------------------------------------
# sbp

mvmr_urates1_sbps2_egfr = MVMR_function("urate_s1","sbp_s2","egfr_sd")
# urate s1 sbp s2
# F-statistic  77.81024  18.80393

mvmr_urates2_sbps1_egfr = MVMR_function("urate_s2","sbp_s1","egfr_sd")
# urate s2 sbp s1
# F-statistic    65.458  17.60337

uratesbp_egfr_mvmr = format_mvmrmeta(mvmr_urates1_sbps2_egfr$result,
                                     mvmr_urates2_sbps1_egfr$result,
                                     "Urate", "SBP")

# format names -----------------------------------------------------------------

uratesbp_egfr_mvmr["exposure"][uratesbp_egfr_mvmr["exposure"] == "Urate (UKB s1) + Urate (UKB s2)"] <- "Urate (UKB Meta)"
uratesbp_egfr_mvmr["exposure"][uratesbp_egfr_mvmr["exposure"] == "SBP (UKB s2) + SBP (UKB s1)"] <- "SBP (UKB Meta)"
uratesbp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

# dbp
mvmr_urates1_dbps2_egfr = MVMR_function("urate_s1","dbp_s2","egfr_sd")
# urate s1 dbp s2
# F-statistic  77.00273  16.79666

mvmr_urates2_dbps1_egfr = MVMR_function("urate_s2","dbp_s1","egfr_sd")
# urate s2 dbp s1
# F-statistic  59.29225  16.55318

uratedbp_egfr_mvmr = format_mvmrmeta(mvmr_urates1_dbps2_egfr$result,mvmr_urates2_dbps1_egfr$result,"Urate","DBP")

# format names -----------------------------------------------------------------

uratedbp_egfr_mvmr["exposure"][uratedbp_egfr_mvmr["exposure"] == "Urate (UKB s1) + Urate (UKB s2)"] <- "Urate (UKB Meta)"
uratedbp_egfr_mvmr["exposure"][uratedbp_egfr_mvmr["exposure"] == "DBP (UKB s2) + DBP (UKB s1)"] <- "DBP (UKB Meta)"
uratedbp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

# pp
mvmr_urates1_pps2_egfr = MVMR_function("urate_s1","pp_s2","egfr_sd")
# urate s1 pp s2
# F-statistic  65.87213  28.22614

mvmr_urates2_pps1_egfr = MVMR_function("urate_s2","pp_s1","egfr_sd")
# urate s2 pp s1
# F-statistic  65.12464  27.35883

uratepp_egfr_mvmr = format_mvmrmeta(mvmr_urates1_pps2_egfr$result,mvmr_urates2_pps1_egfr$result,"Urate", "PP")

# format names -----------------------------------------------------------------

uratepp_egfr_mvmr["exposure"][uratepp_egfr_mvmr["exposure"] == "Urate (UKB s1) + Urate (UKB s2)"] <- "Urate (UKB Meta)"
uratepp_egfr_mvmr["exposure"][uratepp_egfr_mvmr["exposure"] == "PP (UKB s2) + PP (UKB s1)"] <- "PP (UKB Meta)"
uratepp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

# write results ----------------------------------------------------------------

write.table(uratesbp_egfr_mvmr, file = paste0(rdsf_personal,"results/uratesbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(uratedbp_egfr_mvmr, file = paste0(rdsf_personal,"results/uratedbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(uratepp_egfr_mvmr, file = paste0(rdsf_personal,"results/uratepp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
