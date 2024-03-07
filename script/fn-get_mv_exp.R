# Source: https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/

get_mv_exposures <- function(tophits_list, full_gwas_list) {
  
  ###
  ### This is a modified version of `mv_extract_exposures` function in TwoSampleMR package.
  ###
  
  # Collapse list of exposures' tophits into a dataframe
  exposures <- bind_rows(tophits_list)
  
  # clump exposures: this will produce a list of instruments (shared and unique) of the given exposures
  temp <- exposures
  temp$id.exposure <- 1
  temp <- ld_clump_local(temp)
  exposures <- filter(exposures, SNP %in% temp$SNP)
  
  # subset full gwas summary stats of each exposure to the list of SNPs (instruments) produced above
  for (i in 1:length(full_gwas_list)){
    full_gwas_list[[i]] <- full_gwas_list[[i]] %>% filter(SNP %in% exposures$SNP)
  }
  
  # Collapse lists of subset gwas into a dataframe
  d1 <- bind_rows(full_gwas_list) %>%
    distinct()
  
  ###  The logic of next steps is largely unchanged from the original function `mv_extract_exposures`
  
  # get auto-generated ids
  id_exposure <- unique(d1$id.outcome) 
  
  # convert first trait to exposure format  -- exp1 is exposure
  tmp_exposure <- d1 %>% filter(id.outcome == id_exposure[1]) %>% convert_outcome_to_exposure()
  # keep other traits (n>=2) as outcome -- exp2+ are outcomes
  tmp_outcome <- d1 %>% filter(id.outcome != id_exposure[1])
  
  # Harmonise against the first trait
  d <- harmonise_data(exposure_dat = tmp_exposure, 
                      outcome_dat = tmp_outcome, action=2)
  
  # Only keep SNPs that are present in all
  snps_not_in_all <- d %>% 
    dplyr::count(SNP)  %>% 
    filter(n < length(tophits_list)-1) %>%
    pull(SNP)
  d <- filter(d, !SNP %in% snps_not_in_all)
  
  # Subset and concat data
  
  # for exp1 get exposure cols
  dh1x <- d %>% filter(id.outcome == id.outcome[1]) %>% 
    dplyr::select(SNP, dplyr::contains("exposure"))
  # for exp2 get outcome cols
  dh2x <-d %>%  dplyr::select(SNP, dplyr::contains("outcome"))
  # rename outcome to exposure in these
  names(dh2x) <- gsub("outcome", "exposure", names(dh2x) )
  # join together (drop not needed cols)
  exposure_dat <- bind_rows(dh1x, dh2x) %>%  
    dplyr::select(-c("samplesize.exposure" ,"mr_keep.exposure", "pval_origin.exposure")) %>% 
    distinct()
  
  return(exposure_dat)
}

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}

#  function to convert 2SMR format into MVMR format
make_mvmr_input <- function(exposure_dat, outcome.data=""){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  # if (outcome.id.mrbase != "") {
  #   # if mrbase.id is provided
  #   outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
  #                                       outcomes = outcome.id.mrbase)
  # } else if (outcome.data != ""){
  #   # if outcome df is provided
  outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}
