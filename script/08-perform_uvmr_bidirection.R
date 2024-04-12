# Script for performing UVMR analyses

source("code/specify_paths.R")
source("fn-uvmr.R")
source("fn-Isq.R")

# prepare instruments files ----------------------------------------------------

risk_factors <- c("egfr_sd","exurate_sd","urate_clean","sbp_clean","dbp_clean", # for main analyses
                  "urate_clean_s1","urate_clean_s2", # the rest of instruments are for MVMR sensitivity analyses 
                  "sbp_clean_s1","sbp_clean_s2",
                  "dbp_clean_s1","dbp_clean_s2")

all_instruments <- data.frame()
tmp <- NULL

col_order = c("SNP","exposure","effect_allele.exposure","other_allele.exposure",
              "eaf.exposure","beta.exposure","se.exposure","pval.exposure",
              "mr_keep.exposure","pval_origin.exposure","id.exposure",
              "chr.exposure","pos.exposure","samplesize.exposure")

for(i in risk_factors){
  print(i)
  suppressMessages(tmp <- read_tsv(paste0(rdsf_personal,"/data/format_data/",i,"_tophits.tsv")) %>% data.frame())
  tmp <- tmp[,col_order]
  print(nrow(tmp))
  all_instruments = rbind(all_instruments,tmp)
}

egfr2016_exp = extract_instruments("ieu-a-1105",access_token = NULL) %>% dplyr::select(all_of(col_order))
all_instruments = rbind(all_instruments,egfr2016_exp)

all_instruments = F_statistic(all_instruments)

data.table::fwrite(all_instruments, paste0(rdsf_personal,"data/format_data/all_instruments.csv"))

# exposure name is the exposure column in all instruments
# outcome name is the saved files from formatted GWAS file names

# Main analyses between
# urate CKDGen and egfr CKDGen
# urate CKDGen and sbp/dbp UKB 
# sbp/dbp UKB and egfr CKDGen

# Sensitivity analyses between
# urate UKB and eGFR CKDGen

exposure_name <- c("Urate (CKDGen)", "eGFR (CKDGen2019)", "Urate (UKB)", "SBP (UKB)","DBP (UKB)")
outcome_name <- c("egfr_sd", "exurate_sd", "sbp_clean", "dbp_clean","urate_clean")

all_combinations <- expand.grid(exposure = exposure_name, outcome = outcome_name)

filtered_combinations <- subset(all_combinations, !(exposure %in% c("Urate (UKB)","SBP (UKB)","DBP (UKB)") & outcome %in% c("urate_clean","sbp_clean", "dbp_clean") |
                                                    exposure %in% c("urate_clean","sbp_clean", "dbp_clean") & outcome %in% c("Urate (UKB)","SBP (UKB)","DBP (UKB)")))

# define the overlap between strings -------------------------------------------

check_overlap <- function(str1, str2, substrings) {
  for (substr in substrings) {
    if (grepl(substr, str1, ignore.case = TRUE) && grepl(substr, str2, ignore.case = TRUE)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

substrings <- c("sbp", "dbp", "urate", "egfr")

# Create empty results datasets ------------------------------------------------

results_bin <- NULL
plei_bin <- NULL
hetero_bin <- NULL 

# for steiger filtering results ------------------------------------------------

results_sf_bin <- NULL
plei_sf_bin <- NULL
hetero_sf_bin <- NULL

# Do bidirectional MR and save mr, pleio and hetero results --------------------

for (i in 1:nrow(filtered_combinations)) {
  exposure = filtered_combinations[i,"exposure"]
  outcome =  filtered_combinations[i,"outcome"]
  
  print(paste0("Exposure is ", exposure,"; Outcome is ",outcome))
  
    if (!check_overlap(exposure, outcome, substrings)) {
      
      print(paste0("Performing MR for ",exposure, " - ", outcome))
      
      tmp <- uvmr(exposure,outcome)
      results_bin <- rbind(results_bin,tmp[[1]])
      plei_bin <- rbind(plei_bin,tmp[[2]])
      hetero_bin <- rbind(hetero_bin,tmp[[3]])
      results_sf_bin <- rbind(results_sf_bin,tmp[[4]])
      plei_sf_bin <- rbind(plei_sf_bin,tmp[[5]])
      hetero_sf_bin <- rbind(hetero_sf_bin,tmp[[6]])
    }
}

# write results ----------------------------------------------------------------

write.table(results_bin,file = paste0(rdsf_personal,"results/results_bin.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(plei_bin,file = paste0(rdsf_personal,"results/plei_bin.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(hetero_bin,file = paste0(rdsf_personal,"results/hetero_bin.csv"),
            sep= ',', row.names = F,col.names= T)


write.table(results_sf_bin,file = paste0(rdsf_personal,"results/results_sf_bin.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(plei_sf_bin,file = paste0(rdsf_personal,"results/plei_sf_bin.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(hetero_sf_bin,file = paste0(rdsf_personal,"results/hetero_sf_bin.csv"),
            sep= ',', row.names = F,col.names= T)

# additional analyses between bp and egfr 2016----------------------------------
# bp on ieu-a-1105 -------------------------------------------------------------

exposure_name <- c("SBP (UKB)","DBP (UKB)","ieu-a-1105")
outcome_name <- c("sbp_clean", "dbp_clean","ieu-a-1105")

all_combinations <- expand.grid(exposure = exposure_name, outcome = outcome_name)

df1 = uvmr("SBP (UKB)","ieu-a-1105", outcome_sd = 0.24)

df2 = uvmr("DBP (UKB)","ieu-a-1105", outcome_sd = 0.24)

df3 = uvmr("ieu-a-1105", "sbp_clean", exposure_sd = 0.24)

df4 = uvmr("ieu-a-1105","dbp_clean", exposure_sd = 0.24)

sup_results = rbind(df1[[1]],df2[[1]],df3[[1]],df4[[1]])

sup_pleio = rbind(df1[[2]],df2[[2]],df3[[2]],df4[[2]])

sup_hetero = rbind(df1[[3]],df2[[3]],df3[[3]],df4[[3]])

sup_sf_results = rbind(df1[[4]],df2[[4]],df3[[4]],df4[[4]])

sup_sf_pleio = rbind(df1[[5]],df2[[5]],df3[[5]],df4[[5]])

sup_sf_hetero = rbind(df1[[6]],df2[[6]],df3[[6]],df4[[6]])

write.table(sup_results,file = paste0(rdsf_personal,"results/sup_results.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(sup_pleio,file = paste0(rdsf_personal,"results/sup_pleio.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(sup_hetero,file = paste0(rdsf_personal,"results/sup_hetero.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(sup_sf_results,file = paste0(rdsf_personal,"results/sup_sf_results.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(sup_sf_pleio,file = paste0(rdsf_personal,"results/sup_sf_pleio.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(sup_sf_hetero,file = paste0(rdsf_personal,"results/sup_sf_hetero.csv"),
            sep= ',', row.names = F,col.names= T)
