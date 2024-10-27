# Script for UKB GWAS pipeline input for continuous traits

source("fn-ukbgwas-continuous-clean.R")

# read conti data --------------------------------------------------------------
df <- data.table::fread(paste0(rdsf_personal,"data/conti_extract.csv"))
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
# full sample GWAS -------------------------------------------------------------

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

write.table(df_urate_clean, paste0(ukb_input,"df_urate_clean.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_sbp_clean, paste0(ukb_input,"df_sbp_avg_clean.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_dbp_clean, paste0(ukb_input,"df_dbp_avg_clean.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_pp_clean, paste0(ukb_input,"df_pp_avg_clean.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

# Sample split GWAS ------------------------------------------------------------
# Use the complete dataset to split -----------------------------------------------

df_cpt = df[,c("eid","sex","age","sbp_avg","dbp_avg","urate")]
# 502422

df_cpt = merge(df_cpt, link, by = "eid")
# 462826
# remove id thats are identified as outliers in the full sample for all phenotypes
df_cpt = subset(df_cpt, ieu %in% df_urate_clean$IID &
                  ieu %in% df_sbp_clean$IID & 
                  ieu %in% df_dbp_clean$IID &
                  ieu %in% df_pp_clean$IID)
# 413847
sample_id = df_cpt$ieu
# (413847+1)/2 = 206924
set.seed(123)
sample1_id = sample(sample_id, 206924, replace = FALSE)
sample2_id = sample_id[sample_id%notin%sample1_id]

df_urate_clean_s1 = subset(df_urate_clean, IID%in%sample1_id)
# mean 308.7965
# sd 79.78696
# sample size 206924

df_urate_clean_s2 = subset(df_urate_clean, IID%in%sample2_id)
# mean 308.9642
# sd 79.90702
# sample size 206923

df_sbp_avg_clean_s1 = subset(df_sbp_clean, FID%in%df_urate_clean_s1$FID)
# mean 137.744
# sd 18.31552
# sample size 206924

df_sbp_avg_clean_s2 = subset(df_sbp_clean, FID%in%df_urate_clean_s2$FID)
# mean 137.7833
# sd 18.3726
# sample size 206923

df_dbp_avg_clean_s1 = subset(df_dbp_clean, FID%in%df_urate_clean_s1$FID)
# mean 82.12938
# sd 10.04691
# sample size 206924

df_dbp_avg_clean_s2 = subset(df_dbp_clean, FID%in%df_urate_clean_s2$FID)
# mean 82.15788
# sd 10.08448
# sample size 206923

df_pp_avg_clean_s1 = subset(df_pp_clean, FID%in%df_urate_clean_s1$FID)
# mean 55.6146
# sd 13.36189
# sample size 206924

df_pp_avg_clean_s2 = subset(df_pp_clean, FID%in%df_urate_clean_s2$FID)
# mean 55.62538
# sd 13.40611
# sample size 206923

write.table(df_urate_clean_s1, paste0(ukb_input,"df_urate_clean_s1.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)
write.table(df_urate_clean_s2, paste0(ukb_input,"df_urate_clean_s2.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_sbp_avg_clean_s1, paste0(ukb_input,"df_sbp_avg_clean_s1.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)
write.table(df_sbp_avg_clean_s2, paste0(ukb_input,"df_sbp_avg_clean_s2.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_dbp_avg_clean_s1, paste0(ukb_input,"df_dbp_avg_clean_s1.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)
write.table(df_dbp_avg_clean_s2, paste0(ukb_input,"df_dbp_avg_clean_s2.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_pp_avg_clean_s1, paste0(ukb_input,"df_pp_avg_clean_s1.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)
write.table(df_pp_avg_clean_s2, paste0(ukb_input,"df_pp_avg_clean_s2.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)
