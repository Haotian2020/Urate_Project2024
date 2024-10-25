# Script for UKB GWAS pipeline input for continuous traits

source("fn-ukbgwas-continuous-clean.R")

# read conti data --------------------------------------------------------------
df <- data.table::fread("./data/conti_extract.csv")
head(df)
colnames(df) = c("eid","sex","age","sbp1","sbp2","dbp1","dbp2","urate")

df$sbp_avg <- ifelse(
  is.na(df$sbp1) & is.na(df$sbp2), 
  NA, 
  rowMeans(df[, c("sbp1", "sbp2")], na.rm = TRUE)
)

df$dbp_avg <- ifelse(
  is.na(df$dbp1) & is.na(df$dbp2), 
  NA, 
  rowMeans(df[, c("dbp1", "dbp2")], na.rm = TRUE)
)

df$pp_avg = df$sbp_avg - df$dbp_avg

# read linker ------------------------------------------------------------------
link <- data.table::fread("./data/linker_app15825_withexcl.csv")
colnames(link) <- c("ieu","eid")
head(link)

# use continuous data cleaning function to clean each phenotype ----------------

df_urate_clean = singleconti_clean(df,"urate")
# urate mean 308.834600400485
# urate sd 79.895501563483
# urate samplesize 440466

df_sbp_clean = singleconti_clean(df,"sbp_avg")
# sbp mean 137.861762634246
# sbp sd 18.4833994911214
# sbp samplesize 435859

df_dbp_clean = singleconti_clean(df,"dbp_avg")
# dbp mean 82.1555713017935
# dbp sd 10.0886399141079
# dbp samplesize 436083

df_pp_clean = singleconti_clean(df,"pp_avg")
# pp mean 55.6526896881128
# pp sd 13.4170610930985
# pp samplesize 435478

# write gwas input files -------------------------------------------------------

write.table(df_urate_clean,"./data/df_urate_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_sbp_clean,"./data/df_sbp_avg_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_dbp_clean,"./data/df_dbp_avg_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_pp_clean,"./data/df_pp_avg_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

