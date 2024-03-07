# perform uvmr
# source: https://github.com/venexia/T2DMediationMR/blob/master/code/fn-uvmr.R

uvmr <- function(exposure, outcome,ncase = NULL,ncontrol = NULL) {
  
  # Load instruments -----------------------------------------------------------
  
  instruments <- data.table::fread(paste0(rdsf_personal,"data/format_data/all_instruments.csv"), data.table = FALSE)
  
  ## Create exposure dataset ---------------------------------------------------
  
  print(paste0("Reading ",exposure," from all instruments file"))
  tmp_exp <- exposure
  exp <- subset(instruments,exposure == tmp_exp)
  unique(exp$exposure)
  # Extract outcome data -------------------------------------------------------
  
  if (outcome %in% c("egfr_sd","exurate_sd","urate_clean","sbp_clean","dbp_clean","stroke","early50","late60","hpt")) {
    
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
                                             proxies = TRUE,
                                             access_token = NULL)
    
  }
  
  if (!is.null(out)) {
    
    
    ## Harmonise data ----------------------------------------------------------
    
    dat <- suppressMessages(TwoSampleMR::harmonise_data(exposure_dat = exp, 
                                                        outcome_dat = out,action = 2))
    
    ## Perform MR --------------------------------------------------------------
    
    set.seed(123)
    mr <- suppressMessages(TwoSampleMR::mr(dat = dat))
    
    ## Calculate Isq -----------------------------------------------------------
    
    dat_isq <- exp[exp$SNP %in% dat[dat$mr_keep==TRUE,]$SNP,]
    
    isq <- Isq(dat_isq$beta.exposure,dat_isq$se.exposure)
    
    
    
    # Perform Egger intercept test ---------------------------------------------
    
    print("Perform Egger intercept test")
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    plei$Isq <- isq
    
    # Perform heterogeneity test -----------------------------------------------
    
    print("Perform heterogeneity test")
    hetero <- TwoSampleMR::mr_heterogeneity(dat)
    
    # Perform steiger filtering ------------------------------------------------
    
    if(is.null(ncase)&is.null(ncontrol)){
      
      print("The units of both exposure and outcome are in SD")
      dat$units.exposure <- "SD"
      dat$units.outcome <- "SD"
      
      dat <- steiger_filtering(dat)
      print(summary(dat$steiger_dir))
      
    }else
    {
      
      print("The unit of exposure is in SD and the outcome are binary")
      dat$ncase.outcome = ncase
      dat$ncontrol.outcome = ncontrol
      dat$prevalence.outcome = ncase/(ncase + ncontrol)
      
      dat$units.exposure <- "SD"
      dat$units.outcome <- "log odds"

      dat = add_rsq(dat)
      dat$r.outcome = get_r_from_lor(lor = dat$beta.outcome,
                                     af = dat$eaf.exposure,
                                     ncase = unique(dat$ncase.outcome),
                                     ncontrol = unique(dat$ncontrol.outcome),
                                     prevalence = unique(dat$ncase.outcome)/(unique(dat$ncase.outcome+dat$ncontrol.outcome)))
      dat$rsq.outcome = dat$r.outcome^2
      
      dat <- steiger_filtering(dat)
      print(summary(dat$steiger_dir))
    }
    mr_sf = mr(subset(dat, steiger_dir))
    
    mr_hetero_sf <- mr_heterogeneity(subset(dat, steiger_dir))
    
    mr_pleio_sf <- mr_pleiotropy_test(subset(dat, steiger_dir))
  } else {
    
    mr <- NULL
    plei <- NULL
    
  }
  
  r <- list(mr, plei, hetero, mr_sf, mr_hetero_sf, mr_pleio_sf)
  
  return(r)
}
