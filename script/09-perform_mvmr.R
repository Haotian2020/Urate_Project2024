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

mvmr_exurate_bmi_egfr = MVMR_function("exurate_sd","ieu-a-2","egfr_sd")

mvmr_urate_bmi_egfr = MVMR_function("urate","ieu-a-2","egfr_sd")

uvmr("ieu-a-2", "egfr_sd")

# id.exposure id.outcome           outcome        exposure
# 1     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index
# 2     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index
# 3     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index
# 4     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index
# 5     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index
# method nsnp           b         se      pval type
# 1                  MR Egger   79  0.04312029 0.06056701 0.4786513  Ori
# 2           Weighted median   79  0.02144951 0.02240052 0.3382920  Ori
# 3 Inverse variance weighted   79  0.01982116 0.02486343 0.4253342  Ori
# 4               Simple mode   79 -0.03873009 0.05606928 0.4917706  Ori
# 5             Weighted mode   79  0.02385707 0.03135630 0.4490477  Ori
# 
# [[2]]
# id.exposure id.outcome           outcome                      exposure
# 1     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index || id:ieu-a-2
# egger_intercept          se      pval       Isq type
# 1   -0.0007093823 0.001679708 0.6739649 0.9834903  Ori
# 
# [[3]]
# id.exposure id.outcome           outcome                      exposure
# 1     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index || id:ieu-a-2
# 2     ieu-a-2     RyWaLH eGFR (CKDGen2019) Body mass index || id:ieu-a-2
# method        Q Q_df       Q_pval type
# 1                  MR Egger 314.0143   77 2.237977e-30  Ori
# 2 Inverse variance weighted 314.7416   78 3.453153e-30  Ori

MVMR_function("exurate_sd","ieu-a-2","egfr_cys_sd")
# $result
# id.exposure                      exposure id.outcome     outcome nsnp
# 1     ieu-a-2 Body mass index || id:ieu-a-2     DE7oY6 egfr_cys_sd   64
# 2      LrjzHo                Urate (CKDGen)     DE7oY6 egfr_cys_sd   61
# b         se        pval method
# 1 -0.2147843 0.08459193 0.011114950   MVMR
# 2 -0.2370906 0.07336442 0.001230631   MVMR

MVMR_function("urate","ieu-a-2","egfr_cys_sd")

# $result
# id.exposure                      exposure id.outcome     outcome nsnp
# 1     ieu-a-2 Body mass index || id:ieu-a-2     QjQwGj egfr_cys_sd   42
# 2      jQPGtw                   Urate (UKB)     QjQwGj egfr_cys_sd  204
# b         se         pval method
# 1 -0.2404377 0.08047784 2.811535e-03   MVMR
# 2 -0.4076774 0.05746364 1.297933e-12   MVMR

# format names -----------------------------------------------------------------

exurate_sbp_egfr_mvmr = data.frame(mvmr_exurate_sbp_egfr$result)
exurate_sbp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

exurate_dbp_egfr_mvmr = data.frame(mvmr_exurate_dbp_egfr$result)
exurate_dbp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

exurate_pp_egfr_mvmr = data.frame(mvmr_exurate_pp_egfr$result)
exurate_pp_egfr_mvmr$outcome = "eGFR (CKDGen2019)"

# write results ----------------------------------------------------------------
# mvmr results need to be save separately to distinguish the duplicated exposure effects

write.table(exurate_sbp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_sbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(exurate_dbp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_dbp_egfr_mvmr.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(exurate_pp_egfr_mvmr, file = paste0(rdsf_personal,"results/exurate_pp_egfr_mvmr.csv"),
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
