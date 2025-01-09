# function for performing mvmr
# acknowledgment: https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/

getwd()
ao <- available_outcomes()
`%notin%` <- Negate(`%in%`)
source("fn-get_mv_exp.R")
source("fn-ld_clump_local.R")

MVMR_function  <- function(exp1,exp2,outcome1){
  # extract instruments --------------------------------------------------------
  # identify if it is from IEU open GWAS database ------------------------------
  
  if(exp1%in%ao$id & exp2%in%ao$id){
    exptophits1 = extract_instruments(exp1)%>% dplyr::select(c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure",
                                                        "eaf.exposure","beta.exposure","se.exposure","pval_origin.exposure","mr_keep.exposure"))
    exptophits2 = extract_instruments(exp2)%>% dplyr::select(c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure",
                                                        "eaf.exposure","beta.exposure","se.exposure","pval_origin.exposure","mr_keep.exposure"))
    
    tophits_list <- list(exptophits1, exptophits2)
    tophits <- bind_rows(tophits_list) %>% dplyr::pull(SNP)
    
    expgwas1 = extract_outcome_data(snps = tophits, outcomes = exp1, proxies = T) 
    expgwas2 = extract_outcome_data(snps = tophits, outcomes = exp2, proxies = T) 
    
    full_gwas_list <- list(expgwas1, expgwas2)
    print("both expsures are extracted from IEU open GWAS")
    
  }else if(exp1%notin%ao$id & exp2%notin%ao$id){
    
    exptophits1 = fread(paste0(rdsf_personal,"data/format_data/",exp1,"_tophits.tsv"))%>% dplyr::select(c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure",
                                                                                                   "eaf.exposure","beta.exposure","se.exposure","pval.exposure","pval_origin.exposure","mr_keep.exposure"))
    exptophits2 = fread(paste0(rdsf_personal,"data/format_data/",exp2,"_tophits.tsv"))%>% dplyr::select(c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure",
                                                                                                   "eaf.exposure","beta.exposure","se.exposure","pval.exposure","pval_origin.exposure","mr_keep.exposure"))
    
    tophits_list <- list(exptophits1, exptophits2)
    tophits <- bind_rows(tophits_list) %>% dplyr::pull(SNP)
    
    expgwas1 = vroom(paste0(rdsf_personal,"data/format_data/",exp1,"_GWAS_tidy_outcome.csv"))
    expgwas2 = vroom(paste0(rdsf_personal,"data/format_data/",exp2,"_GWAS_tidy_outcome.csv"))
    full_gwas_list <- list(expgwas1, expgwas2)
    
    print("both expsures are extracted from local")
    
  }else if(exp1%notin%ao$id & exp2%in%ao$id){
    
    exptophits1 = fread(paste0(rdsf_personal,"data/format_data/",exp1,"_tophits.tsv")) %>% dplyr::select(c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure",
                                                                                                    "eaf.exposure","beta.exposure","se.exposure","pval.exposure","pval_origin.exposure","mr_keep.exposure"))
    exptophits2 = extract_instruments(exp2) %>% dplyr::select(c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure",
                                                         "eaf.exposure","beta.exposure","se.exposure","pval.exposure","pval_origin.exposure","mr_keep.exposure"))

    exptophits2 = exptophits2[,colnames(exptophits1)]
    
    tophits_list <- list(exptophits1, exptophits2)
    tophits <- bind_rows(tophits_list) %>% dplyr::pull(SNP)
    
    expgwas1 = vroom(paste0(rdsf_personal,"data/format_data/",exp1,"_GWAS_tidy_outcome.csv")) %>% dplyr::select(-c("chr.outcome", "pos.outcome", "pval_origin.outcome"))
    expgwas2 = extract_outcome_data(snps = tophits, outcomes = exp2, proxies = T) %>% dplyr::select(colnames(expgwas1))
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
    
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome1)
    
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
  res_bmis <- data.frame(res_bmis$result)
  res_bmis$method = 'MVMR'
  print(res_bmis)
  
  mvmr_input_dat <- MVMR::format_mvmr(
    BXGs = mvdat$exposure_beta,
    seBXGs = mvdat$exposure_se,
    BYG = mvdat$outcome_beta,
    seBYG = mvdat$outcome_se
  )
  
  sres <- strength_mvmr(r_input = mvmr_input_dat, gencov = 0)
  print(sres)
  
  return(list(result = res_bmis, F_stat = sres))
}

