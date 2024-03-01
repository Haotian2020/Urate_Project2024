source("fn-ukbgwas-binary-output-format.R")

# read and format output from gwas files ---------------------------------------

hpt = BinaryOutputFormatFunction("df_hpt_imputed.txt.gz","Hypertension",133680,329146)

early50 = BinaryOutputFormatFunction("df_early_imputed.txt.gz","Early-onset hypertension (50)",6934,329146)

late60 = BinaryOutputFormatFunction("df_late_imputed.txt.gz","Late-onset hypertension (60)",95583,329146)

# write formatted files for binary output --------------------------------------

write.table(hpt, file = paste0(rdsf_personal,"data/format_data/hpt_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(early50, file = paste0(rdsf_personal,"data/format_data/early50_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(late60, file = paste0(rdsf_personal,"data/format_data/late60_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)
