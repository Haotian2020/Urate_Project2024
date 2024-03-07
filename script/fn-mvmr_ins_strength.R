mvmr_instruments_strength <- function(exp1,exp2,out1){
  
  # Load instruments -----------------------------------------------------------
  
  tophits1 = vroom(paste0(rdsf_personal,"data/format_data/",exp1,"_tophits.tsv"))
  tophits2 = vroom(paste0(rdsf_personal,"data/format_data/",exp2,"_tophits.tsv"))
  
  tophits_list <- list(tophits1, 
                       tophits2)
  
  tophits <- bind_rows(tophits_list) %>% pull(SNP)

  # Load two exposure gwas -----------------------------------------------------
  expgwas1 = vroom(paste0(rdsf_personal,"data/format_data/",exp1,"_GWAS_tidy_outcome.csv"))
  expgwas2 = vroom(paste0(rdsf_personal,"data/format_data/",exp2,"_GWAS_tidy_outcome.csv"))
  
  full_gwas_list <- list(expgwas1, expgwas2)
  
  exposure_dat <- get_mv_exposures(tophits_list, full_gwas_list)

  outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
                                   filename =  paste0(rdsf_personal,"data/format_data/",out1,"_GWAS_tidy_outcome.csv"),
                                   sep = ',',
                                   snp_col = "SNP",
                                   beta_col = "beta.outcome",
                                   se_col = "se.outcome",
                                   eaf_col = "eaf.outcome",
                                   effect_allele_col = "effect_allele.outcome",
                                   other_allele_col = "other_allele.outcome",
                                   pval_col = "pval.outcome")
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
  mvmr_input_dat <- MVMR::format_mvmr(BXGs = mvdat$exposure_beta,
                                seBXGs = mvdat$exposure_se,
                                BYG = mvdat$outcome_beta,
                                seBYG = mvdat$outcome_se)
  sres <- strength_mvmr(r_input = mvmr_input_dat, gencov = 0)
  return(sres)
}

