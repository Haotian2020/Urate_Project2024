source("specify_paths.R")
getwd()
`%notin%` <- Negate(`%in%`)

df <- data.table::fread(paste0(rdsf_path,"linker_app15825.csv"),
                        header = TRUE, 
                        data.table = FALSE)

eur <- data.table::fread(paste0(rdsf_personal,"data/GWAS_input/data.europeans.qctools.txt"),
                         header = FALSE, 
                         data.table = FALSE)

df <- df[(df$ieu %in% eur$V1),]

# Recommended exclusions ------------------------------------------------------
# for sex_mismatch, putative_sex_chromosome_aneuploidy and het_missing_outliers

excl <- data.table::fread(paste0(rdsf_personal,"data/GWAS_input/data.combined_recommended.qctools.txt"),
                          header = FALSE, 
                          data.table = FALSE)

df <- df[!(df$ieu %in% excl$V1),]

# Remove relateds -------------------------------------------------------------
rel_high <- data.table::fread(paste0(rdsf_personal,"data/GWAS_input/data.highly_relateds.qctools.txt"),
                              header = FALSE, 
                              data.table = FALSE)

df <- df[!(df$ieu %in% rel_high$V1),]

# Remove withdrals
ukb_withdraw <- data.table::fread(paste0(ukb_withdraw,"/w15825_2023-04-25.csv"),
                                  header = FALSE,
                                  data.table = FALSE)

df <- df[!(df$ieu %in% ukb_withdraw$V1),]

# Save final linker file -------------------------------------------------------
data.table::fwrite(df,paste0(rdsf_personal,"data/linker_app15825_withexcl.csv"))
