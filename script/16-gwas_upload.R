# website source: https://mrcieu.github.io/GwasDataImport/articles/import_pipeline_new.html

# install.packages("remotes")
# remotes::install_github("MRCIEU/GwasDataImport")

library(GwasDataImport)

# make sure you have a token to access OpenGWAS
# readRenviron(".Renviron")

upload_gwas <- function(file_path, metadata, columns) {
  x <- GwasDataImport::Dataset$new(filename = file_path)
  x$collect_metadata(metadata)
  x$api_metadata_upload()
  print(paste("GWAS ID:", x$igd_id))
  print(paste("Filename:", x$filename))
  
  cat("Pausing before determining columns. Perform your necessary external operations.\n")
  readline(prompt = "Press Enter to continue once ready...")
  
  x$determine_columns(columns)
  x$format_dataset()
  x$metadata_uploaded <- TRUE
  x$api_gwasdata_upload()
  x$delete_wd()
}

# Urate GWAS
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/urate_imputed.txt.gz"),
  metadata = list(
    trait = "Urate",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Risk factor",
    subcategory = "Biomarker",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    sample_size = 440466,
    author = "Tang H",
    year = 2024,
    unit = "SD"
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# Systolic Blood Pressure (SBP) GWAS
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/sbp_avg_imputed.txt.gz"),
  metadata = list(
    trait = "Systolic Blood Pressure (SBP)",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Risk factor",
    subcategory = "Blood pressure",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    sample_size = 435859,
    author = "Tang H",
    year = 2024,
    unit = "SD"
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# Diastolic Blood Pressure (DBP) GWAS
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/dbp_avg_imputed.txt.gz"),
  metadata = list(
    trait = "Diastolic Blood Pressure (DBP)",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Risk factor",
    subcategory = "Blood pressure",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    sample_size = 436083,
    author = "Tang H",
    year = 2024,
    unit = "SD"
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# Pulse Pressure (PP) GWAS
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/pp_avg_imputed.txt.gz"),
  metadata = list(
    trait = "Pulse Pressure (PP)",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Risk factor",
    subcategory = "Blood pressure",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    sample_size = 435478,
    author = "Tang H",
    year = 2024,
    unit = "SD"
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# ns ukb gwas
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/df_ns_imputed.txt.gz"),
  metadata = list(
    trait = "Nephrotic syndrome",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Binary",
    subcategory = "Kidney",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    ncase = 349,
    ncontrol = 459687,
    author = "Tang H",
    year = 2024,
    unit = "beta_bolt",
    sample_size = 460036
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# hypertension ukb gwas
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/df_hpt_imputed.txt.gz"),
  metadata = list(
    trait = "Hypertension",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Binary",
    subcategory = "Blood pressure",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    ncase = 133680,
    ncontrol = 329146,
    author = "Tang H",
    year = 2024,
    unit = "beta_bolt",
    sample_size = 462826
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# early onset
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/df_early_imputed.txt.gz"),
  metadata = list(
    trait = "Early-onset hypertension",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Binary",
    subcategory = "Blood pressure",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    ncase = 6934,
    ncontrol = 329146,
    author = "Tang H",
    year = 2024,
    unit = "beta_bolt",
    sample_size = 336080
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)

# late onset
upload_gwas(
  file_path = paste0(rdsf_personal,"data/GWAS_imputed/df_late_imputed.txt.gz"),
  metadata = list(
    trait = "Late-onset hypertension",
    group_name = "public",
    build = "HG19/GRCh37",
    category = "Binary",
    subcategory = "Blood pressure",
    ontology = "NA",
    population = "European",
    sex = "Males and Females",
    ncase = 95583,
    ncontrol = 329146,
    author = "Tang H",
    year = 2024,
    unit = "beta_bolt",
    sample_size = 424729
  ),
  columns = list(
    chr_col = "CHR",
    pos_col = "BP",
    ea_col = "ALLELE1",
    oa_col = "ALLELE0",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P_BOLT_LMM_INF",
    snp_col = "SNP",
    eaf_col = "A1FREQ"
  )
)
