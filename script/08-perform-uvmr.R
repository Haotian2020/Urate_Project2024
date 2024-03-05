source("code/specify_paths.R")
source("fn-uvmr.R")
source("fn-steiger.R")
source("fn-Isq.R")

# prepare instruments files ----------------------------------------------------

risk_factors <- c("egfr_sd","exurate_sd","urate_clean","sbp_clean","dbp_clean")

all_instruments <- data.frame()
tmp <- NULL

col_order = c("SNP","exposure","effect_allele.exposure","other_allele.exposure",
              "eaf.exposure","beta.exposure","se.exposure","pval.exposure",
              "mr_keep.exposure","pval_origin.exposure","id.exposure",
              "chr.exposure","pos.exposure","samplesize.exposure")

for(i in risk_factors){
  print(i)
  tmp <- read_tsv(paste0(rdsf_personal,"/data/format_data/",i,"_tophits.tsv")) %>% data.frame()
  tmp <- tmp[,col_order]
  print(head(tmp))

  all_instruments = rbind(all_instruments,tmp)
}

all_instruments = F_statistic(all_instruments)

data.table::fwrite(all_instruments, paste0(rdsf_personal,"data/format_data/all_instruments.csv"))

# Create empty results datasets ------------------------------------------------

results_bin <- NULL #data.table::fread("output/results.csv", data.table = FALSE)
plei_bin <- NULL #data.table::fread("output/plei.csv", data.table = FALSE)
hetero_bin <- NULL 

results_sf_bin <- NULL
plei_sf_bin <- NULL
hetero_sf_bin <- NULL

# exposure name is the exposure column in all instruments
# outcome name is the saved files from formatted GWAS file names

# Main analyses between 
# urate CKDGen and egfr CKDGen
# urate CKDGen and sbp/dbp UKB 
# sbp/dbp UKB and egfr CKDGen

# Sensitivity analyses between
# urate UKB and eGFR CKDGen

exposure_name <- c("Urate CKDGen", "eGFR CKDGen", "urate_clean", "sbp_clean","dbp_clean")
outcome_name <- c("egfr_sd", "exurate_sd", "sbp_clean", "dbp_clean","urate_clean")

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

# Do bidirectional MR and save mr, pleio and hetero results --------------------

for (exposure in exposure_name) {
  for (outcome in outcome_name) {
    if (!check_overlap(exposure, outcome, substrings)) {
      print(paste0("Performing MR for ",exposure, "-", outcome))
      tmp <- uvmr(exposure,outcome)
      results_bin <- rbind(results_bin,tmp[[1]])
      plei_bin <- rbind(plei_bin,tmp[[2]])
      hetero_bin <- rbind(hetero_bin,tmp[[3]])
      results_sf_bin <- rbind(results_sf_bin,tmp[[4]])
      plei_sf_bin <- rbind(plei_sf_bin,tmp[[5]])
      hetero_sf_bin <- rbind(hetero_sf_bin,tmp[[6]])
    }
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
