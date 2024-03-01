# continuous data cleaning function --------------------------------------------

singleconti_clean <- function(dat,var){
  print(paste0("total data set No.row is ",nrow(dat)))
  dat = select(dat,c("eid",all_of(var))) %>% na.omit()
  print(paste0("after removing NA, No.row is ",nrow(dat)))
  print(colnames(dat))
  dat = merge(dat, link, by = "eid")
  print(paste0("after merging with link, No.row is ",nrow(dat)))
  dat = dat[,c(3,3,2)]
  print(head(dat))
  colnames(dat) = c("FID","IID",var)
  var_mean = mean(dat[[var]])
  var_sd = sd(dat[[var]])
  # print(paste0("mean of ",var," is ",var_mean))
  # print(paste0("sd of ",var," is ",var_sd ))
  dat_clean = subset(dat,dat[[3]]<=(var_mean+4*var_sd) & dat[[3]]>=(var_mean - 4*var_sd))
  print(paste0("the number of outliers were removed ", nrow(dat) - nrow(dat_clean)))
  print(paste0("mean of ",var," is ",mean(dat_clean[[3]])))
  print(paste0("sd of ",var," is ",sd(dat_clean[[3]])))
  print(paste0("cleaned data No.row is ",nrow(dat_clean)))
  return(dat_clean)
}

# continuous gwas data format function -----------------------------------------
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