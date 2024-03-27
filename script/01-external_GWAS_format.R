rm(list=ls())
graphics.off()
`%notin%` <- Negate(`%in%`)

# Load libraries ---------------------------------------------------------------

source("lib setting.R")

# Specify paths ----------------------------------------------------------------

source("specify_paths.R")

# Format external Urate GWAS ---------------------------------------------------

exurate_gwas = vroom(paste0(rdsf_personal,"data/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt")) %>% 
  subset(!is.na(RSID))

exurate_gwas_format = format_data(
  exurate_gwas,
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
  samplesize_col = "n_total_sum"
)
exurate_gwas_format$outcome = "Urate (CKDGen)"

# Transform to SD unit ---------------------------------------------------------

exurate_sd_gwas_format = exurate_gwas_format
exurate_sd_gwas_format$beta.outcome = exurate_sd_gwas_format$beta.outcome/1.5
exurate_sd_gwas_format$se.outcome = exurate_sd_gwas_format$se.outcome/1.5

# Find tophits -----------------------------------------------------------------

exurate_sd_tophits = ld_clump_local(exurate_sd_gwas_format)

# Save formated outcome and instruments ----------------------------------------

write_tsv(exurate_sd_tophits, file = paste0(rdsf_personal, 'data/format_data/exurate_sd_tophits.tsv'))
write.table(exurate_sd_gwas_format, file = paste0(rdsf_personal,"data/format_data/exurate_sd_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)


# Format eGFR_crea 2019 GWAS ---------------------------------------------------
egfr_gwas_outcome = fread(paste0(rdsf_personal,"data/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt")) %>% data.frame()

egfr_gwas_outcome_format = format_data(egfr_gwas_outcome,
                                       type = 'outcome',
                                       snp_col = "RSID",
                                       beta_col = "Effect",
                                       se_col = "StdErr",
                                       effect_allele_col = "Allele1",
                                       other_allele_col = "Allele2",
                                       eaf_col = "Freq1",
                                       pval_col = "P.value",
                                       samplesize_col = "n_total_sum",
                                       chr_col = "Chr",
                                       pos_col = "Pos_b37")

egfr_gwas_outcome_format$outcome = "eGFR (CKDGen2019)"

# Transform to SD unit ---------------------------------------------------------
egfr_gwas_outcome_sd_format = egfr_gwas_outcome_format
egfr_gwas_outcome_sd_format$beta.outcome = egfr_gwas_outcome_sd_format$beta.outcome/0.13
egfr_gwas_outcome_sd_format$se.outcome = egfr_gwas_outcome_sd_format$se.outcome/0.13

# Find tophits -----------------------------------------------------------------
egfr_sd_tophits = ld_clump_local(egfr_gwas_outcome_sd_format)

# Save formated outcome and instruments ----------------------------------------
write_tsv(egfr_sd_tophits, file = paste0(rdsf_personal,'data/format_data/egfr_sd_tophits.tsv'))
write.table(egfr_gwas_outcome_sd_format, file = paste0(rdsf_personal,'data/format_data/egfr_sd_GWAS_tidy_outcome.csv'),
            sep= ',', row.names = F,col.names= T)


# Format stroke GWAS -----------------------------------------------------------

stroke_gwas_outcome = vroom(paste0(rdsf_personal,"data/MEGASTROKE.1.AS.EUR.out"))
stroke_gwas_outcome$pval = stroke_gwas_outcome$`P-value`

stroke_gwas_outcome_format = format_data(stroke_gwas_outcome, type = "outcome",
                                         snp_col = "MarkerName",
                                         beta_col = "Effect",
                                         se_col = "StdErr",
                                         effect_allele_col = "Allele1",
                                         other_allele_col = "Allele2",
                                         eaf_col = "Freq1",
                                         pval_col = "pval")
stroke_gwas_outcome_format$outcome =  'Stroke'

stroke_gwas_outcome_format$ncontrol.outcome = 406111
stroke_gwas_outcome_format$ncase.outcome = 40585

# Save formated outcome --------------------------------------------------------
write.table(stroke_gwas_outcome_format, file = paste0(rdsf_personal,'data/format_data/stroke_GWAS_tidy_outcome.csv'),
            sep= ',', row.names = F,col.names= T)
