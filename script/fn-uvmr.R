# perform uvmr
# source: https://github.com/venexia/T2DMediationMR/blob/master/code/fn-uvmr.R

uvmr <- function(exposure, outcome) {
  
  # Load instruments -----------------------------------------------------------
  
  instruments <- data.table::fread(paste0(rdsf_personal,"data/format_data/all_instruments.csv"), data.table = FALSE)
  
  ## Create exposure dataset ---------------------------------------------------
  
  print(paste0("Reading ",exposure," from all instruments file"))
  tmp_exp <- exposure
  exp <- subset(instruments,exposure == tmp_exp)

  # Extract outcome data -------------------------------------------------------
  
  if (outcome %in% c("egfr_sd","exurate_sd","urate_clean","sbp_clean","dbp_clean")) {
    
    print(paste0("Reading ",outcome," from local formatted data"))
    out <- TwoSampleMR::read_outcome_data(snps = exp$SNP,
                             filename = paste0(rdsf_personal,"data/format_data/",outcome,"_GWAS_tidy_outcome.csv"),
                             sep = ",",
                             phenotype_col = "outcome",
                             snp_col = "SNP",
                             beta_col = "beta.outcome",
                             se_col = "se.outcome",
                             eaf_col = "eaf.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             pval_col = "pval.outcome",
                             samplesize_col = "samplesize.outcome")
  } else {
    
    print("Reading from IEU open GWAS database")
    out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP,
                                             outcomes = outcome,
                                             proxies = FALSE,
                                             access_token = NULL)
    
  }
  
  if (!is.null(out)) {
    
    
    ## Harmonise data ----------------------------------------------------------
    
    dat <- suppressMessages(TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                                        outcome_dat = out,action = 2))
    
    ## Perform MR --------------------------------------------------------------
    
    results <- suppressMessages(TwoSampleMR::mr(dat = dat))
    
    ## Calculate Isq -----------------------------------------------------------
    
    dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
    
    isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
    
    results$Isq <- isq
    
    # Perform Egger intercept test ---------------------------------------------
    
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    # Perform heterogeneity test -----------------------------------------------
    
    hetero <- TwoSampleMR::mr_heterogeneity(dat)
    
  } else {
    
    results <- NULL
    plei <- NULL
    
  }
  
  r <- list(results, plei,hetero)
  
  return(r)
}