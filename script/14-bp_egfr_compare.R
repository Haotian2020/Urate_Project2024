# use bp instruments to compare the SNP effects in two egfr gwas 
# Wuttke, Li Y., Li M., Sieber, Feitosa, Gorski et al., 2019: eGFR, BUN and CKD associations; trans-ethnic and European American ancestry
# Gorski et al., 2015: GWAS of kidney function decline in individuals of European descent

for(i in c("sbp_clean","dbp_clean")){
  
  exp = data.frame(read_tsv(paste0(rdsf_personal,"data/format_data/",i,"_tophits.tsv")))
  
  egfr_out1 <- read_outcome_data(
    snps = exp$SNP,
    filename = paste0(rdsf_personal,"data/format_data/egfr_GWAS_tidy_outcome.csv"),
    sep = ",",
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    samplesize_col = "samplesize.outcome"
  ) %>% convert_outcome_to_exposure()
  
  egfr_out1$exposure = "eGFR Wuttke et al 2019"
  
  egfr_out2 = extract_outcome_data(snps = egfr_out1$SNP,outcomes = "ieu-a-1105",access_token = NULL)
  
  
}

