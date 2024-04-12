# function for formatting UKB GWAS summary level data for continuous trait

ContiOutputFormatFunction = function(dat_name,outcome_name,mean,sd,samplesize){
  print("reading data from rdsf")
  dat <-vroom(paste0(rdsf_personal,"data/GWAS_imputed/",dat_name))
  print(colnames(dat))
  if("P_BOLT_LMM"%in%colnames(dat)){
    print("using P_BOLT_LMM as p value")
    format_dat = format_data(dat, type = "outcome",
                             snp_col = "SNP",
                             beta_col = "BETA",
                             se_col = "SE",
                             effect_allele_col = "ALLELE1",
                             other_allele_col = "ALLELE0",
                             eaf_col = "A1FREQ",
                             pval_col = "P_BOLT_LMM",
                             chr_col = "CHR",
                             pos_col = "BP")
  }else{
    print("using P_BOLT_LMM_INF as p value")
    format_dat = format_data(dat, type = "outcome",
                             snp_col = "SNP",
                             beta_col = "BETA",
                             se_col = "SE",
                             effect_allele_col = "ALLELE1",
                             other_allele_col = "ALLELE0",
                             eaf_col = "A1FREQ",
                             pval_col = "P_BOLT_LMM_INF",
                             chr_col = "CHR",
                             pos_col = "BP")
    
  }
  print("adding sample size and outcome name")
  format_dat$outcome = outcome_name
  format_dat$samplesize.outcome = samplesize
  print(colnames(format_dat))
  print("standardize the unit (sd unit)")
  format_dat$beta.outcome = format_dat$beta.outcome/sd
  format_dat$se.outcome = format_dat$se.outcome/sd
  print(head(format_dat))
  return(format_dat)
}