# Script for performing UVMR analyses

source("code/specify_paths.R")
source("fn-uvmr.R")
source("fn-Isq.R")

# prepare instruments files ----------------------------------------------------

risk_factors <- c("egfr_sd","exurate_sd","urate","sbp","dbp", "pp", # for main analyses
                  "urate_s1","urate_s2", # the rest of instruments are for MVMR sensitivity analyses 
                  "sbp_s1","sbp_s2",
                  "dbp_s1","dbp_s2",
                  "pp_s1","pp_s2")

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

egfr2016_exp = extract_instruments("ieu-a-1105") %>% dplyr::select(all_of(col_order))
egfr2016_exp$exposure = "eGFR (CKDGen2016)"

all_instruments = rbind(all_instruments, egfr2016_exp)

all_instruments = F_statistic(all_instruments)

data.table::fwrite(all_instruments, paste0(rdsf_personal,"data/format_data/all_instruments.csv"))

# exposure name is the exposure column in all instruments
# outcome name is the saved files from formatted GWAS file names

# Main analyses between
# urate CKDGen and egfr CKDGen
# urate CKDGen and sbp/dbp/pp UKB 
# sbp/dbp/pp UKB and egfr CKDGen

# Sensitivity analyses between
# urate UKB and eGFR CKDGen

exposure_name <- c("exurate_sd", "egfr_sd", "urate", "sbp","dbp", "pp")
outcome_name <- c("egfr_sd", "exurate_sd", "sbp", "dbp", "pp", "urate")

all_combinations <- expand.grid(exposure = exposure_name, outcome = outcome_name)

filtered_combinations <- subset(all_combinations, !(exposure %in% c("urate","sbp","dbp", "pp") & outcome %in% c("urate","sbp", "dbp", "pp") |
                                                    exposure %in% c("urate","sbp", "dbp", "pp") & outcome %in% c("urate","sbp", "dbp", "pp")))

# define the overlap between strings -------------------------------------------

check_overlap <- function(str1, str2, substrings) {
  for (substr in substrings) {
    if (grepl(substr, str1, ignore.case = TRUE) && grepl(substr, str2, ignore.case = TRUE)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

substrings <- c("sbp", "dbp", "pp", "urate", "egfr")

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
df1 = uvmr("sbp", "ieu-a-1105", outcome_sd = 0.24)

df2 = uvmr("dbp", "ieu-a-1105", outcome_sd = 0.24)

df3 = uvmr("pp", "ieu-a-1105", exposure_sd = 0.24)

df4 = uvmr("ieu-a-1105", "sbp", exposure_sd = 0.24)

df5 = uvmr("ieu-a-1105","dbp", exposure_sd = 0.24)

df6 = uvmr("ieu-a-1105","pp", exposure_sd = 0.24)

sup_results <- do.call(rbind, lapply(list(df1, df2, df3, df4, df5, df6), `[[`, 1))
sup_pleio <- do.call(rbind, lapply(list(df1, df2, df3, df4, df5, df6), `[[`, 2))
sup_hetero <- do.call(rbind, lapply(list(df1, df2, df3, df4, df5, df6), `[[`, 3))
sup_sf_results <- do.call(rbind, lapply(list(df1, df2, df3, df4, df5, df6), `[[`, 4))
sup_sf_pleio <- do.call(rbind, lapply(list(df1, df2, df3, df4, df5, df6), `[[`, 5))
sup_sf_hetero <- do.call(rbind, lapply(list(df1, df2, df3, df4, df5, df6), `[[`, 6))

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
