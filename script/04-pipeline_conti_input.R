source("fn-ukbgwas-continuous-clean.R")

# read conti data --------------------------------------------------------------
df <- data.table::fread("./data/conti_extract.csv")
head(df)
colnames(df) = c("eid","sex","age","sbp","dbp","urate")

# read linker ------------------------------------------------------------------
link <- data.table::fread("./data/linker_app15825_withexcl.csv")
colnames(link) <- c("ieu","eid")
head(link)

# use continuous data cleaning function to clean each phenotype ----------------

df_urate_clean = singleconti_clean(df,"urate")
# urate mean 308.834600400485
# urate sd 79.895501563483
# urate samplesize 440466

df_sbp_clean = singleconti_clean(df,"sbp")
# sbp mean 139.78564403065
# sbp sd 19.5241937805439
# sbp samplesize 432099

df_dbp_clean = singleconti_clean(df,"dbp")
# dbp mean 82.1392934230647
# dbp sd 10.6322476240864
# dbp samplesize 432253

# write gwas input files -------------------------------------------------------

write.table(df_urate_clean,"./data/df_urate_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_sbp_clean,"./data/df_sbp_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_dbp_clean,"./data/df_dbp_clean.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)
