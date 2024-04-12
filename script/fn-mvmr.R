# function for performing mvmr

getwd()
ao <- available_outcomes(access_token = NULL)
`%notin%` <- Negate(`%in%`)
source("fn-get_mv_exp.R")
source("fn-ld_clump_local.R")

MVMR_function  <- function(exp1,exp2,outcome1){
  
  # extract instruments --------------------------------------------------------
  # identify if it is from IEU open GWAS database ------------------------------
  
  if(exp1%in%ao$id & exp2%in%ao$id){
    exptophits1 = extract_instruments(exp1, access_token = NULL)%>% select(-c("chr.exposure","pos.exposure"))
    exptophits2 = extract_instruments(exp2, access_token = NULL)%>% select(-c("chr.exposure","pos.exposure"))
    
    tophits_list <- list(exptophits1, exptophits2)
    tophits <- bind_rows(tophits_list) %>% pull(SNP)
    
    expgwas1 = extract_outcome_data(snps = tophits, outcomes = exp1, proxies = T, access_token = NULL) 
    expgwas2 = extract_outcome_data(snps = tophits, outcomes = exp2, proxies = T, access_token = NULL) 
    
    full_gwas_list <- list(expgwas1, expgwas2)
    print("both expsures are extracted from IEU open GWAS")
    
  }else if(exp1%notin%ao$id & exp2%notin%ao$id){
    
    exptophits1 = read_tsv(paste0(rdsf_personal,"data/format_data/",exp1,"_tophits.tsv"))%>% select(-c("chr.exposure","pos.exposure"))
    exptophits2 = read_tsv(paste0(rdsf_personal,"data/format_data/",exp2,"_tophits.tsv"))%>% select(-c("chr.exposure","pos.exposure"))
    
    tophits_list <- list(exptophits1, exptophits2)
    tophits <- bind_rows(tophits_list) %>% pull(SNP)
    
    expgwas1 = vroom(paste0(rdsf_personal,"data/format_data/",exp1,"_GWAS_tidy_outcome.csv"))
    expgwas2 = vroom(paste0(rdsf_personal,"data/format_data/",exp2,"_GWAS_tidy_outcome.csv"))
    full_gwas_list <- list(expgwas1, expgwas2)
    
    print("both expsures are extracted from local")
    
  }else if(exp1%notin%ao$id & exp2%in%ao$id){
    
    exptophits1 = read_tsv(paste0(rdsf_personal,"data/format_data/",exp1,"_tophits.tsv")) %>% select(-c("chr.exposure","pos.exposure"))
    exptophits2 = extract_instruments(exp2, access_token = NULL) %>% select(-c("chr.exposure","pos.exposure"))

    print(colnames(exptophits1))
    print(colnames(exptophits1))
    
    exptophits2 = exptophits2[,colnames(exptophits1)]
    
    tophits_list <- list(exptophits1, exptophits2)
    tophits <- bind_rows(tophits_list) %>% pull(SNP)
    
    expgwas1 = vroom(paste0(rdsf_personal,"data/format_data/",exp1,"_GWAS_tidy_outcome.csv")) %>% select(-c("chr.outcome", "pos.outcome", "pval_origin.outcome"))
    expgwas2 = extract_outcome_data(snps = tophits, outcomes = exp2,access_token = NULL, proxies = T) %>% select(colnames(expgwas1))
    full_gwas_list <- list(expgwas1, expgwas2)
    print("exp1 is extracted from local")
    print("exp2 is extracted from IEU open GWAS")
  }
  
  # make instruments list ------------------------------------------------------
  
  print("create the whole exp infor from two exposures")
  exposures <- bind_rows(tophits_list)
  print("using get_mv_exp function")
  exposure_dat <- get_mv_exposures(tophits_list, full_gwas_list)
  print("exposure data prepared")
  
  # identify if the outcome is from IEU open GWAS database ---------------------
  
  if(outcome1%in%ao$id){
    
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome1,access_token = NULL)
    
  }else{
    
    outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
                                     filename =  paste0(rdsf_personal,"data/format_data/",outcome1,"_GWAS_tidy_outcome.csv"),
                                     sep = ",",
                                     snp_col = "SNP",
                                     beta_col = "beta.outcome",
                                     se_col = "se.outcome",
                                     eaf_col = "eaf.outcome",
                                     effect_allele_col = "effect_allele.outcome",
                                     other_allele_col = "other_allele.outcome",
                                     pval_col = "pval.outcome")%>%mutate(outcome = outcome1)
  }
  
  # perform mvmr ---------------------------------------------------------------
  
  mvdat <- TwoSampleMR::mv_harmonise_data(exposure_dat, outcome_dat)
  res_bmis <- TwoSampleMR::mv_multiple(mvdat)
  return(res_bmis)
}


