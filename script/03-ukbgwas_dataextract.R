# Scripts for extracting data from UKB, including ID, sex, age, genotyping chip,
# BP, urate and ICD10

# read and write urate, sbp and dbp data from ukb ------------------------------

vars <- c("eid","31-0.0","21022-0.0","4080-0.0","4079-0.0", "6177-0.0",# ID, sex, age, SBP, DBP,
          "30880-0.0") # urate

df <- fread(paste0(rdsf_path,"data.48733.csv"),
            sep = ",", 
            header = TRUE, 
            select = vars,
            data.table = FALSE)

data.table::fwrite(df,paste0(rdsf_personal,"data/conti_extract.csv"))

# read and write birth infor and icd10 data from ukb ---------------------------

vars <- c("eid","34-0.0","52-0.0")# year of birth, month of birth

df_birth <- fread(paste0(rdsf_path,"data.48733.csv"),
                  sep = ",", 
                  header = TRUE, 
                  select = vars,
                  data.table = FALSE)

vars <- c(paste0("41270-0.",seq(0,225)), # ICD10 diagnoses
          paste0("41280-0.",seq(0,225))) # date of ICD10 diagnoses

df_icd10 <- fread(paste0(rdsf_path,"data.48733.csv"),
                  sep = ",", 
                  header = TRUE, 
                  select = vars,
                  data.table = FALSE)

df_icd10time <- cbind(df_birth,df_icd10)

data.table::fwrite(df_icd10time,paste0(rsdf_personal,"data/GWAS_input/icd10withtime_extract.csv"))

# create covariates file ---------------------------------------------------------

covs = fread(paste0(rdsf_personal,"data/data.covariates.bolt.txt"))

vars = c("eid","31-0.0","21022-0.0")
df_age <- fread(paste0(rdsf_path,"data.48733.csv"),
                sep = ",", 
                header = TRUE, 
                select = vars,
                data.table = FALSE)
colnames(df_age) = c("eid","sex","age")

link <- data.table::fread(paste0(rdsf_personal,"./data/linker_app15825_withexcl.csv"))
colnames(link) <- c("ieu","eid")
df_age = merge(df_age,link,by="eid")
df_age = df_age[,c("ieu","sex","age")]
colnames(df_age) = c("IID","sex","age")

covs_age = merge(covs,df_age,by = "IID")
covs_age = covs_age[,c("FID","IID", "sex.x","chip","age")]
colnames(covs_age) = c("FID","IID", "sex","chip","age")

# write covariates (sex, chip, age) file ----------------------------------------

data.table::fwrite(covs_age,paste0(rsdf_personal,"data/data.covariates_age.bolt.txt"))
