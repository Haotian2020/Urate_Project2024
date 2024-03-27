# use bp instruments to compare the SNP effects in two egfr gwas 
# Wuttke, Li Y., Li M., Sieber, Feitosa, Gorski et al., 2019: eGFR, BUN and CKD associations; trans-ethnic and European American ancestry
# Pattaro, Teumer, Gorski, Chu, Li, and Mijatovic et al., 2016: Genetic associations at 53 loci highlight cell types and biological pathways relevant for kidney function

for(i in c("SBP (UKB)","DBP (UKB)")){
  
  exp = data.frame(fread(paste0(rdsf_personal,"data/format_data/all_instruments.csv"))) %>% subset(exposure == i)
  
  egfr_out1 <- read_outcome_data(
    snps = exp$SNP,
    filename = paste0(rdsf_personal,"data/format_data/egfr_GWAS_tidy_outcome.csv"),
    sep = ",",
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    samplesize_col = "samplesize.outcome"
  ) %>% convert_outcome_to_exposure()
  
  egfr_out1$exposure = "log(eGFR)"
  egfr_out2 = extract_outcome_data(snps = egfr_out1$SNP,outcomes = "ieu-a-1105",proxies = F, access_token = NULL)
  egfr_out2$outcome = "log(eGFR)"
  # egfr_out2$beta.outcome = egfr_out2$beta.outcome/0.24
  # egfr_out2$se.outcome = egfr_out2$se.outcome/0.24
  egfr_har = harmonise_data(exposure_dat = egfr_out1,outcome_dat = egfr_out2,action = 2)
  
  df <- egfr_har[,c("SNP","beta.exposure","se.exposure","exposure","beta.outcome","se.outcome","outcome")]
  print(nrow(df))

  p <- ggplot2::ggplot(data=df, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(ggplot2::aes(intercept=0, slope=1), show.legend=TRUE) +
    ggplot2::xlim(min(df$beta.exposure - df$se.exposure,df$beta.outcome - df$se.outcome),max(df$beta.outcome + df$se.outcome,df$beta.exposure + df$se.exposure)) + 
    ggplot2::ylim(min(min(df$beta.exposure - df$se.exposure,df$beta.outcome - df$se.outcome)),max(df$beta.outcome + df$se.outcome,df$beta.exposure + df$se.exposure)) +
    ggplot2::scale_colour_manual(values=c("black", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
    ggplot2::labs(x = paste0("SNP effect on ", unique(egfr_har$exposure), " (CKDGen2019)"),
                  y = paste0("SNP effect on ", unique(egfr_har$outcome), " (CKDGen2016)"))  +
    ggplot2::theme(legend.position="top", legend.direction="vertical") +
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))

  pdf(paste0(rdsf_personal,"results/",i," on 2 egfr scatter.pdf"),width = 6, height = 6)
  print(p)
  dev.off()
}


# calculate the SD of log(eGFR) for Pattaro 2016 -------------------------------

egfr2016cohort = read_excel(path = paste0(rdsf_personal,"data/1011.xlsx"),sheet =1)
egfr2016cohort = egfr2016cohort[2:nrow(egfr2016cohort),c("Study", "Sample size","Age Mean (SD)" ,"eGFRcrea Mean(SD) ml / min /\r\n1.73m2")]
colnames(egfr2016cohort) = c("study","samplesize","age","egfr")
egfr_split <- strsplit(gsub("\\)", "", gsub("\\(", " ", egfr2016cohort$egfr)), " ")
egfr_mean <- sapply(egfr_split, function(x) as.numeric(x[1]))
egfr_sd <- sapply(egfr_split, function(x) as.numeric(x[2]))

egfr2016cohort$egfr_mean <- as.numeric(egfr_mean)
egfr2016cohort$egfr_sd <- as.numeric(egfr_sd)
egfr2016cohort$samplesize <- as.numeric(egfr2016cohort$samplesize)
egfr2016cohort <- egfr2016cohort[egfr2016cohort$study != "PREVEND", ]

egfr2016cohort$log_egfr = log(egfr2016cohort$egfr_mean)

egfr2016cohort$log_egfr_sd = sqrt(log(egfr2016cohort$egfr_sd^2 / egfr2016cohort$egfr_mean^2 + 1))

total_log_sd <- sqrt(sum(egfr2016cohort$log_egfr_sd^2 * egfr2016cohort$samplesize) / sum(egfr2016cohort$samplesize))

# total_log_sd = 0.2407372
