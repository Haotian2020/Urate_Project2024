method_order <- factor(
  levels = c(
    "Inverse variance weighted",
    "Weighted median",
    "Steiger Filtering",
    "Weighted mode",
    "MR Egger",
    "Simple mode")
)

risk_factor_order <- factor(
  levels = c(
    "Urate (CKDGen)",
    "eGFR (CKDGen)",
    "Urate (UKB)",
    "SBP (UKB)",
    "DBP (UKB)",
    "Urate (UKB s1)", "Urate (UKB s2)", "SBP (UKB s1)",  
    "SBP (UKB s2)", "DBP (UKB s1)", "DBP (UKB s2)")
)

apply_custom_order <- function(df, order_type) {
  if (order_type == "risk_factor_order") {
    
    ordered_df <- df[order(factor(df$exposure, levels = levels(risk_factor_order))), ]
    
  } else if (order_type == "method_order") {

    ordered_df <- df[order(factor(df$method, levels = levels(method_order))), ]
    
  } else if(order_type == "both"){
    
    ordered_df <- df[order(factor(df$exposure, levels = levels(risk_factor_order)),
                                     factor(df$method, levels = levels(method_order))), ]
  }
  
  return(ordered_df)
}

# genetic instruments sheet ----------------------------------------------------

dat = fread(paste0(rdsf_personal,"data/format_data/all_instruments.csv"))

names(dat) <- sub("\\.exposure", "", names(dat))

selected_col = c("exposure","SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize","F_stat","chr","pos")

dat <- dat %>% select(all_of(selected_col))

dat <- apply_custom_order(dat,"risk_factor_order")

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee1.csv"))

# hetero test sheet ------------------------------------------------------------

dat1 = fread(paste0(rdsf_personal,"results/hetero_bin.csv"))

dat2 = fread(paste0(rdsf_personal,"results/hpt_hetero.csv"))

dat3 = fread(paste0(rdsf_personal,"results/pos_hetero.csv"))

dat = rbind(dat1,dat2,dat3) %>% select(-c("id.exposure","id.outcome")) %>% 
  filter(type == "Ori") %>% split_outcome() %>% select(c("exposure", "outcome",	"method",	"Q",	"Q_df",	"Q_pval"))

dat <- apply_custom_order(dat,order_type = "risk_factor_order")

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee2.csv"))

# pleiotropy sheet -------------------------------------------------------------

dat1 = fread(paste0(rdsf_personal,"results/plei_bin.csv"))

dat2 = fread(paste0(rdsf_personal,"results/hpt_pleio.csv"))

dat3 = fread(paste0(rdsf_personal,"results/pos_pleio.csv"))

dat = rbind(dat1,dat2,dat3) %>% select(-c("id.exposure","id.outcome")) %>% filter(type == "Ori")

dat <- dat[,c("exposure","outcome","egger_intercept","se","pval","Isq")]

dat <- apply_custom_order(dat,"risk_factor_order")

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee3.csv"))

# bidirectional MR sheet -------------------------------------------------------

dat = fread(paste0(rdsf_personal,"results/results_bin.csv")) %>% select(-c("type","id.exposure","id.outcome"))

dat <- dat[,c("exposure","outcome","method",	"nsnp",	"b", "se", "pval")]

dat <- apply_custom_order(dat,"risk_factor_order")

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee4.csv"))

# bidirectional MR sheet (Steiger) ---------------------------------------------

dat = fread(paste0(rdsf_personal,"results/results_sf_bin.csv")) %>% select(-c("type","id.exposure","id.outcome"))

dat <- dat[,c("exposure","outcome","method",	"nsnp",	"b", "se", "pval")]

dat <- apply_custom_order(dat,"risk_factor_order")

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee5.csv"))

# binary outcomes --------------------------------------------------------------

dat = fread(paste0(rdsf_personal,"results/hpt_res.csv")) %>% 
  select(-c("id.exposure","id.outcome")) %>% 
  filter(type == "Ori") %>%
  select(c("exposure","outcome","method",	"nsnp",	"b", "se", "pval","type")) %>% split_outcome()

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee6_1.csv"))

dat = fread(paste0(rdsf_personal,"results/hpt_res.csv")) %>% 
  select(-c("id.exposure","id.outcome")) %>% 
  filter(type == "Steiger") %>%
  select(c("exposure","outcome","method",	"nsnp",	"b", "se", "pval","type")) %>% split_outcome()

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee6_2.csv"))

dat = fread(paste0(rdsf_personal,"results/pos_res.csv")) %>% 
  select(-c("id.exposure","id.outcome")) %>% 
  filter(type == "Ori") %>%
  select(c("exposure","outcome","method",	"nsnp",	"b", "se", "pval","type")) %>% split_outcome()

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee6_3.csv"))

dat = fread(paste0(rdsf_personal,"results/pos_res.csv")) %>% 
  select(-c("id.exposure","id.outcome")) %>% 
  filter(type == "Steiger") %>%
  select(c("exposure","outcome","method",	"nsnp",	"b", "se", "pval","type")) %>% split_outcome()

data.table::fwrite(dat,paste0(rdsf_personal,"results/shee6_4.csv"))

# mvmr -------------------------------------------------------------------------

# directly use the resutls from "09-perform_mvmr.R"





