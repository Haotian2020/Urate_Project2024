# function for formatting UKB GWAS summary level data for binary trait

BinaryOutputFormatFunction = function(dat_name,outcome_name,Ncase,Ncontrol){
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
  format_dat$ncase.outcome = Ncase
  format_dat$ncontrol.outomce = Ncontrol
  format_dat$samplesize.outcome = Ncase + Ncontrol
  print(colnames(format_dat))
  mu = Ncase/(Ncase + Ncontrol)
  print(paste0("prevalence of ",outcome_name," is ", mu))
  print("calculate the lambda")
  x = mu*(1-mu)
  print("standardize the unit to mu")
  format_dat$beta.outcome = format_dat$beta.outcome/x
  format_dat$se.outcome = format_dat$se.outcome/x
  print(head(format_dat))
  print("No tophits being extracted")
  return(format_dat)
}