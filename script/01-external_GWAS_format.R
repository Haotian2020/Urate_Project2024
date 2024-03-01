rm(list=ls())
graphics.off()

# Load libraries ---------------------------------------------------------------
source("lib setting.R")

# Specify paths ----------------------------------------------------------------
source("specify_paths.R")

# Format external Urate GWAS ---------------------------------------------------
exurate_gwas = vroom(paste0(rdsf_personal,"data/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt")) %>% 
  subset(!is.na(RSID))

# Format data ------------------------------------------------------------------
exurate_gwas_format = format_data(
  exurate_gwas_keep,
  type = "outcome",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
  chr_col = "Chr",
  pos_col = "Pos_b37",
  samplesize_col = "n_total_sum")
exurate_gwas_format$outcome = "Urate CKDGen"

# Transform to SD unit ---------------------------------------------------------
exurate_sd_gwas_format = exurate_gwas_format
exurate_sd_gwas_format$beta.outcome = exurate_sd_gwas_format$beta.outcome/1.5
exurate_sd_gwas_format$se.outcome = exurate_sd_gwas_format$se.outcome/1.5

# Find tophits -----------------------------------------------------------------
exurate_sd_tophits = 
  exurate_gwas_format %>% 
  filter(pval.outcome < 5e-8) %>% 
  convert_outcome_to_exposure() %>% 
  clump_data(., clump_r2 = 0.001)

exurate_samplesizeinfo = subset(exurate_gwas_format,SNP%in%exurate_sd_tophits$SNP)[,c("SNP","samplesize.outcome")]
colnames(exurate_samplesizeinfo) = c("SNP","samplesize.exposure")
exurate_sd_tophits = merge(exurate_sd_tophits,exurate_samplesizeinfo,by = c("SNP"))
# Save formated outcome and instruments -----------------------------------------------------
write_tsv(exurate_sd_tophits, path = paste0(rdsf_personal, 'data/format_data/exurate_sd_tophits.tsv'))
write.table(exurate_sd_gwas_format, file = paste0(rdsf_personal,"data/format_data/exurate_sd_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

# Format eGFR crea GWAS --------------------------------------------------------



















