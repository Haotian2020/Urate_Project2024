# Script for formatting GWAS summary level from UKB pipeline output of binary traits

source("fn-ukbgwas-binary-output-format.R")

# Read and format output from gwas files ---------------------------------------

hpt = BinaryOutputFormatFunction("df_hpt_imputed.txt.gz","Hypertension (UKB)",133680,329146)

early50 = BinaryOutputFormatFunction("df_early_imputed.txt.gz","Early-onset hypertension (UKB)",6934,329146)

late60 = BinaryOutputFormatFunction("df_late_imputed.txt.gz","Late-onset hypertension (UKB)",95583,329146)

# Write formatted files for binary output --------------------------------------

data.table::fwrite(hpt,paste0(rdsf_personal,"data/format_data/hpt_GWAS_tidy_outcome.csv"))

data.table::fwrite(early50,paste0(rdsf_personal,"data/format_data/early50_GWAS_tidy_outcome.csv"))

data.table::fwrite(late60,paste0(rdsf_personal,"data/format_data/late60_GWAS_tidy_outcome.csv"))
