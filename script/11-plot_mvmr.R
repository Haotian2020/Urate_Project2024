# Script for plotting MVMR figures

source("format_meta.R")

# urate ckdgen + bp ukb on egfr ckdgen -----------------------------------------
# compare with urate ckdgen on egfr ckdgen and bp ukb on egfr ckdgen -----------

# Load data --------------------------------------------------------------------

exurate_sbp_egfr_mvmr = fread(paste0(rdsf_personal,"results/exurate_sbp_egfr_mvmr.csv"))

exurate_dbp_egfr_mvmr = fread(paste0(rdsf_personal,"results/exurate_dbp_egfr_mvmr.csv"))

exurate_pp_egfr_mvmr = fread(paste0(rdsf_personal,"results/exurate_pp_egfr_mvmr.csv"))

res = fread(paste0(rdsf_personal,"results/results_bin.csv"))

exposures <- c("SBP", "DBP", "PP")
labels <- c("A", "B", "C")

# Iterate through the exposures, labels, and file names
for (i in seq_along(exposures)) {
  exp <- exposures[i]
  label <- labels[i]

  # Combine data for the current exposure
  mydata <- rbind(
    res %>% subset(exposure == "Urate (CKDGen)" & outcome == "eGFR (CKDGen2019)" & method == "Inverse variance weighted"),
    res %>% subset(exposure == paste0(exp, " (UKB)") & outcome == "eGFR (CKDGen2019)" & method == "Inverse variance weighted"),
    get(paste0("exurate_", tolower(exp), "_egfr_mvmr")),
    fill = TRUE
  )

  # Generate odds ratios and format the data
  mydata <- mydata %>% generate_odds_ratios() %>%
    mutate(method = if_else(method == "Inverse variance weighted", "IVW", method))

  mydata$beta <- format(round(mydata$b, digits = 3), nsmall = 3)
  mydata$CI <- paste0(format(round(mydata$lo_ci, digits = 3), nsmall = 3), ", ",
                      format(round(mydata$up_ci, digits = 3), nsmall = 3))
  mydata$pvalue <- format(round(mydata$pval, digits = 3))
  mydata$method <- factor(mydata$method, levels = c("IVW", "MVMR"))

  mydata$exposure <- factor(
    mydata$exposure,
    levels = c("Urate (CKDGen)", "SBP (UKB)", "DBP (UKB)", "PP (UKB)", "eGFR(CKDGen)")
  )

  sorted_index <- order(mydata$exposure, mydata$method)
  mydata <- mydata[sorted_index, ]

  if (nrow(mydata) >= 4) {
    for (j in c(seq(2, nrow(mydata), 4), seq(3, nrow(mydata), 4), seq(4, nrow(mydata), 4))) {
      mydata$outcome[j] <- NA
    }
  }

  if (nrow(mydata) >= 2) {
    for (j in seq(2, nrow(mydata), 2)) {
      mydata$exposure[j] <- NA
    }
  }

  mydata <- data.frame(mydata)

  # Prepare table text for the forest plot
  tabletext <- cbind(
    c("Exposure", as.character(mydata[, 'exposure'])),
    c("Outcome", as.character(mydata[, 'outcome'])),
    c("Approach", as.character(mydata[, 'method'])),
    c('Number of SNPs', as.character(mydata[, 'nsnp'])),
    c("Beta", as.character(mydata[, 'beta'])),
    c("95% CI", as.character(mydata[, 'CI'])),
    c("p-value", as.character(mydata[, 'pvalue']))
  )

  # Generate the forest plot
  p <- forestplot(
    tabletext,
    graph.pos = 4,
    mean = as.numeric(rbind(NA, cbind(mydata[, 'beta']))),
    lower = as.numeric(rbind(NA, cbind(mydata[, 'lo_ci']))),
    upper = as.numeric(rbind(NA, cbind(mydata[, 'up_ci']))),
    new_page = FALSE,
    psignif = 0.05,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    hrzl_lines = list("4" = gpar(lty = 1, lwd = 1, col = "black")),
    boxsize = 0.3,
    line.margin = 0.2,
    lty.ci = 1,
    col = fpColors(box = "black", lines = "darkgray"),
    lwd.ci = 2,
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    graphwidth = unit(150, "mm"),
    is.summary = c(TRUE, rep(FALSE, nrow(tabletext))),
    colgap = unit(5, "mm"),
    zero = 0.0,
    xticks = c(-0.2, -0.1, 0, 0.1),
    clip = c(-0.2, 0.1),
    xlab = paste0("Beta (with 95% CI) for SD-unit of log(eGFR) per SD unit change in exposure")
  )

  # Save the plot to a PDF
  output_file <- paste0(
  "results/exurate ", tolower(exp), " ukb on egfr forestplot.pdf")
  pdf_file <- paste0(rdsf_personal, output_file)
  pdf(pdf_file, width = 18, height = 2.4)
  plot.new()
  print(p)
  mtext(paste0(label, ")"), side = 3, line = 2, adj = 0, cex = 1.5, padj = 0)
  dev.off()
}

# urate ukb  + bp ukb (sample split method) --> egfr ckdgen --------------------
# compare with uvmr results, including urate ukb on egfr ckdgen and bp ukb on egfr ckdgen
# here we need to compare with the results from sample-split method ------------
# doing uvmr first and then meta

urate_egfr_mr_meta = format_meta(uvmr("urate_s1", "egfr_sd")[[1]],
                                 uvmr("urate_s2", "egfr_sd")[[1]])

sbp_egfr_mr_meta = format_meta(uvmr("sbp_s1", "egfr_sd")[[1]],
                               uvmr("sbp_s2", "egfr_sd")[[1]])

dbp_egfr_mr_meta = format_meta(uvmr("dbp_s1", "egfr_sd")[[1]],
                               uvmr("dbp_s2", "egfr_sd")[[1]])

pp_egfr_mr_meta = format_meta(uvmr("pp_s1", "egfr_sd")[[1]],
                               uvmr("pp_s2", "egfr_sd")[[1]])

# format results ---------------------------------------------------------------

df <- rbind(sbp_egfr_mr_meta, dbp_egfr_mr_meta, urate_egfr_mr_meta, pp_egfr_mr_meta)%>% subset(method == "Inverse variance weighted")

df$outcome <- "eGFR (CKDGen2019)"
df$exposure <- c("SBP (UKB Meta)", "DBP (UKB Meta)", "Urate (UKB Meta)", "PP (UKB Meta)")

# Define the exposures and their corresponding output labels
exposures <- c("SBP", "DBP", "PP")
labels <- c("A", "B", "C")

# Iterate through the exposures
for (i in seq_along(exposures)) {
  exp <- exposures[i]
  label <- labels[i]

  # Load mvmr data ---------------------------------------------------------------
  mvmr_file <- fread(paste0(rdsf_personal, "results/urate", tolower(exp), "_egfr_mvmr.csv"))

  # Filter and bind data dynamically
  mydata <- rbind(
    df %>% subset(exposure == "Urate (UKB Meta)"),
    df %>% subset(exposure == paste0(exp, " (UKB Meta)")),
    mvmr_file
  )

  # Generate odds ratios and format data
  mydata <- mydata %>% generate_odds_ratios() %>%
    mutate(method = if_else(method == "Inverse variance weighted", "IVW", method))

  mydata$beta <- format(round(mydata$b, digits = 3), nsmall = 3)
  mydata$CI <- paste0(format(round(mydata$lo_ci, digits = 3), nsmall = 3), ", ",
                      format(round(mydata$up_ci, digits = 3), nsmall = 3))
  mydata$pvalue <- format(round(mydata$pval, digits = 3))
  mydata$method <- factor(
    mydata$method, 
    levels = c("IVW", "MVMR"))
  mydata$exposure <- factor(
    mydata$exposure, 
    levels = c("Urate (UKB Meta)", 
    "SBP (UKB Meta)", 
    "DBP (UKB Meta)", 
    "PP (UKB Meta)", 
    "eGFR (CKDGen2019)"))

  sorted_index <- order(mydata$exposure, mydata$method)
  mydata <- mydata[sorted_index, ]
  mydata <- data.frame(mydata)

  for (j in c(seq(2, nrow(mydata), 4),
            seq(3, nrow(mydata), 4),
            seq(4, nrow(mydata), 4))) {
  mydata$outcome[j] <- NA
}

for (j in c(seq(2, nrow(mydata), 2))) {
  mydata$exposure[j] <- NA
}

  mydata <- data.frame(mydata)

  # Prepare table text for the forest plot
  tabletext <- cbind(
    c("Exposure", as.character(mydata[, 'exposure'])),
    c("Outcome", as.character(mydata[, 'outcome'])),
    c("Approach", as.character(mydata[, 'method'])),
    c('Number of SNPs', as.character(mydata[, 'nsnp'])),
    c("Beta", as.character(mydata[, 'beta'])),
    c("95% CI", as.character(mydata[, 'CI'])),
    c("p-value", as.character(mydata[, 'pvalue']))
  )

  # Generate the forest plot
  p <- forestplot(
    tabletext,
    graph.pos = 4,
    mean = as.numeric(rbind(NA, cbind(mydata[, 'beta']))),
    lower = as.numeric(rbind(NA, cbind(mydata[, 'lo_ci']))),
    upper = as.numeric(rbind(NA, cbind(mydata[, 'up_ci']))),
    new_page = FALSE,
    psignif = 0.05,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    hrzl_lines = list("4" = gpar(lty = 1, lwd = 1, col = "black")),
    boxsize = 0.3,
    line.margin = 0.2,
    lty.ci = 1,
    col = fpColors(box = "black", lines = "darkgray"),
    lwd.ci = 2,
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    graphwidth = unit(150, "mm"),
    is.summary = c(TRUE, rep(FALSE, nrow(tabletext))),
    colgap = unit(5, "mm"),
    zero = 0.0,
    xticks = c(-0.2, -0.1, 0, 0.1),
    clip = c(-0.2, 0.1),
    xlab = paste0("Beta (with 95% CI) for SD-unit of log(eGFR) per SD unit change in exposure")
  )

  # Save the plot to a PDF
  pdf_file <- paste0(rdsf_personal, "results/urate meta ", tolower(exp), " ukb on egfr forestplot.pdf")
  pdf(pdf_file, width = 18, height = 2.4)
  plot.new()
  print(p)
  mtext(paste0(label, ")"), side = 3, line = 2, adj = 0, cex = 1.5, padj = 0)
  dev.off()
}
