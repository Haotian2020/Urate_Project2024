# Script for plotting MVMR figures

source("format_meta.R")

# urate ckdgen + bp ukb on egfr ckdgen -----------------------------------------
# compare with urate ckdgen on egfr ckdgen and bp ukb on egfr ckdgen -----------

# Load data --------------------------------------------------------------------

exurate_sbp_egfr_mvmr = fread(paste0(rdsf_personal,"results/exurate_sbp_egfr_mvmr.csv"))

exurate_dbp_egfr_mvmr = fread(paste0(rdsf_personal,"results/exurate_dbp_egfr_mvmr.csv"))

res = fread(paste0(paste0(rdsf_personal,"results/results_bin.csv")))

# format results ---------------------------------------------------------------

replacement_rules <- c("urate_clean" = "Urate (UKB)", "sbp_clean" = "SBP (UKB)", "dbp_clean" = "DBP (UKB)", "Urate CKDGen" = "Urate (CKDGen)",
                       "eGFR (CKDGen)" = "eGFR (CKDGen2019)")

for (pattern in names(replacement_rules)) {
  res$exposure <- sub(pattern, replacement_rules[pattern], res$exposure)
  res$outcome <- sub(pattern, replacement_rules[pattern], res$outcome)
}

# prepare dataset for figures --------------------------------------------------

mydata = rbind(res %>% subset(exposure == "Urate (CKDGen)" & outcome == "eGFR (CKDGen2019)") %>% subset(method == "Inverse variance weighted"),
               res %>% subset(exposure == "SBP (UKB)" & outcome == "eGFR (CKDGen2019)") %>% subset(method == "Inverse variance weighted"),
               exurate_sbp_egfr_mvmr,fill = T)

mydata = rbind(res %>% subset(exposure == "Urate (CKDGen)" & outcome == "eGFR (CKDGen2019)") %>% subset(method == "Inverse variance weighted"),
               res %>% subset(exposure == "DBP (UKB)" & outcome == "eGFR (CKDGen2019)") %>% subset(method == "Inverse variance weighted"),
               exurate_dbp_egfr_mvmr,fill = T)

mydata = mydata %>% generate_odds_ratios()

mydata$beta <-
  format(round(mydata$b, digits = 3),
         nsmall = 3)
mydata$CI <-
  paste0(format(round(mydata$lo_ci, digits = 3),
                nsmall = 3),
         ", ",
         format(round(mydata$up_ci, digits = 3),
                nsmall = 3))
mydata$pvalue <-
  format(round(mydata$pval, digits = 3))
mydata$method <-
  factor(
    mydata$method,
    levels = c(
      "Inverse variance weighted","MVMR"))
mydata
unique(mydata$exposure)
mydata$exposure <-
  factor(
    mydata$exposure,
    levels = c(
      "Urate (CKDGen)","SBP (UKB)","DBP (UKB)","eGFR(CKDGen)"))

sorted_index <- order(mydata$exposure,mydata$method)
mydata = mydata[sorted_index,]

for (i in c(seq(2, nrow(mydata), 4),
            seq(3, nrow(mydata), 4),
            seq(4, nrow(mydata), 4))) {
  mydata$outcome[i] <- NA
}

for (i in c(seq(2, nrow(mydata), 2))) {
  mydata$exposure[i] <- NA
}

mydata = data.frame(mydata)
mydata
tabletext <- cbind(
  c("Expsoure", as.character(mydata[, 'exposure'])),
  c("Outcome", as.character(mydata[, 'outcome'])),
  c("Approach", as.character(mydata[, 'method'])),
  c('Number of SNPs', as.character(mydata[, 'nsnp'])),
  c("Beta", as.character(mydata[, 'beta'])),
  c("95% CI", as.character(mydata[, 'CI'])),
  c("p-value", as.character(mydata[, 'pvalue'])))
tabletext
max(mydata$up_ci)
min(mydata$lo_ci)

dev.off()
p <- forestplot(
  tabletext,
  graph.pos = 4,
  mean = as.numeric(rbind(NA, cbind(mydata[, 'beta']))),
  lower = as.numeric(rbind(NA, cbind(mydata[, 'lo_ci']))),
  upper = as.numeric(rbind(NA, cbind(mydata[, 'up_ci']))),
  new_page = F,
  psignif = 0.05,
  txt_gp = fpTxtGp(
    label = gpar(cex = 1),
    ticks = gpar(cex = 1),
    xlab = gpar(cex = 1),
    title = gpar(cex = 1)
  ),
  hrzl_lines = list("4" = gpar(lty = 1, lwd = 1, col = "black")),
  
  boxsize = 0.15,
  line.margin = 0.1,
  lty.ci = 1,
  col = fpColors(box = "black", lines = "darkgray"),
  lwd.ci = 1,
  ci.vertices = T,
  ci.vertices.height = 0.15,
  graphwidth = unit(150, "mm"),
  is.summary = c(T, rep(F, nrow(tabletext))),
  colgap = unit (5, "mm"),
  zero = 0.0,
  xticks = c(-0.2,-0.2,-0.1,0,0.1),
  clip = c(-0.2,0.1),
  xlab = paste0("Beta (with 95% CI) for SD-unit of log(eGFR)  per SD unit change in exposure"))
p

# c(-0.3,-0.2,-0.1,0,0.1,0.2)
dev.off()
tiff(paste0(rdsf_personal,"results/exurate sbp ukb on egfr forestplot.tiff"),width = 18, height = 2.5, res = 300, units = "in")
plot.new()
print(p)
mtext("A)",side = 3,line = 2,adj = 0,cex = 1.5,padj = 0)
dev.off()

tiff(paste0(rdsf_personal,"results/exurate dbp ukb on egfr forestplot.tiff"),width = 18, height = 2.5, res = 300, units = "in")
plot.new()
print(p)
mtext("B)",side = 3,line = 2,adj = 0,cex = 1.5,padj = 0)
dev.off()

# urate ukb  + bp ukb (sample split method) --> egfr ckdgen --------------------
# compare with urate ukb on egfr ckdgen and bp ukb on egfr ckdgen --------------

# Load data --------------------------------------------------------------------

uratesbp_egfr_mvmr = fread(paste0(rdsf_personal,"results/uratesbp_egfr_mvmr.csv"))

uratedbp_egfr_mvmr = fread(paste0(rdsf_personal,"results/uratedbp_egfr_mvmr.csv"))

# here we need to compare with the results from sample-split method ------------

urate_egfr_mr_meta = format_meta(uvmr("Urate (UKB s1)","egfr_sd")[[1]],uvmr("Urate (UKB s2)","egfr_sd")[[1]])

sbp_egfr_mr_meta = format_meta(uvmr("SBP (UKB s1)","egfr_sd")[[1]],uvmr("SBP (UKB s2)","egfr_sd")[[1]])

dbp_egfr_mr_meta = format_meta(uvmr("DBP (UKB s1)","egfr_sd")[[1]],uvmr("DBP (UKB s2)","egfr_sd")[[1]])

# format results ---------------------------------------------------------------

df <- rbind(sbp_egfr_mr_meta, dbp_egfr_mr_meta, urate_egfr_mr_meta)%>% subset(method == "Inverse variance weighted")

df$outcome = "eGFR (CKDGen2019)"
df$exposure[1] = "SBP (UKB Meta)"
df$exposure[2] = "DBP (UKB Meta)"
df$exposure[3] = "Urate (UKB Meta)"


mydata = rbind(df %>% subset(exposure == "Urate (UKB Meta)"),
               df %>% subset(exposure == "SBP (UKB Meta)"),
               uratesbp_egfr_mvmr)

mydata = rbind(df %>% subset(exposure == "Urate (UKB Meta)"),
               df %>% subset(exposure == "DBP (UKB Meta)"),
               uratedbp_egfr_mvmr)

mydata = mydata %>% generate_odds_ratios()

mydata$beta <-
  format(round(mydata$b, digits = 3),
         nsmall = 3)
mydata$CI <-
  paste0(format(round(mydata$lo_ci, digits = 3),
                nsmall = 3),
         ", ",
         format(round(mydata$up_ci, digits = 3),
                nsmall = 3))
mydata$pvalue <-
  format(round(mydata$pval, digits = 3))
mydata$method <-
  factor(
    mydata$method,
    levels = c(
      "Inverse variance weighted","MVMR"))
mydata$exposure <-
  factor(
    mydata$exposure,
    levels = c(
      "Urate (UKB Meta)","SBP (UKB Meta)","DBP (UKB Meta)","eGFR (CKDGen2019)"))

sorted_index <- order(mydata$exposure,mydata$method)
mydata = mydata[sorted_index,]

for (i in c(seq(2, nrow(mydata), 4),
            seq(3, nrow(mydata), 4),
            seq(4, nrow(mydata), 4))) {
  mydata$outcome[i] <- NA
}

for (i in c(seq(2, nrow(mydata), 2))) {
  mydata$exposure[i] <- NA
}

mydata = data.frame(mydata)
mydata
tabletext <- cbind(
  c("Expsoure", as.character(mydata[, 'exposure'])),
  c("Outcome", as.character(mydata[, 'outcome'])),
  c("Approach", as.character(mydata[, 'method'])),
  c('Number of SNPs', as.character(mydata[, 'nsnp'])),
  c("Beta", as.character(mydata[, 'beta'])),
  c("95% CI", as.character(mydata[, 'CI'])),
  c("p-value", as.character(mydata[, 'pvalue'])))
tabletext
max(mydata$up_ci)
min(mydata$lo_ci)

dev.off()
p <- forestplot(
  tabletext,
  graph.pos = 4,
  mean = as.numeric(rbind(NA, cbind(mydata[, 'beta']))),
  lower = as.numeric(rbind(NA, cbind(mydata[, 'lo_ci']))),
  upper = as.numeric(rbind(NA, cbind(mydata[, 'up_ci']))),
  new_page = F,
  psignif = 0.05,
  txt_gp = fpTxtGp(
    label = gpar(cex = 1),
    ticks = gpar(cex = 1),
    xlab = gpar(cex = 1),
    title = gpar(cex = 1)
  ),
  hrzl_lines = list("4" = gpar(lty = 1, lwd = 1, col = "black")),
  
  boxsize = 0.15,
  line.margin = 0.1,
  lty.ci = 1,
  col = fpColors(box = "black", lines = "darkgray"),
  lwd.ci = 1,
  ci.vertices = T,
  ci.vertices.height = 0.15,
  graphwidth = unit(150, "mm"),
  is.summary = c(T, rep(F, nrow(tabletext))),
  colgap = unit (5, "mm"),
  zero = 0.0,
  xticks = c(-0.2,-0.2,-0.1,0,0.1),
  clip = c(-0.2,0.1),
  xlab = paste0("Beta (with 95% CI) for SD-unit of log(eGFR)  per SD unit change in exposure"))
p

# c(-0.3,-0.2,-0.1,0,0.1,0.2)
dev.off()
tiff(paste0(rdsf_personal,"results/urate meta sbp ukb on egfr forestplot.tiff"),width = 18, height = 2.5, res = 300, units = "in")
plot.new()
print(p)
mtext("A)",side = 3,line = 2,adj = 0,cex = 1.5,padj = -1)
dev.off()

tiff(paste0(rdsf_personal,"results/urate meta dbp ukb on egfr forestplot.tiff"),width = 18, height = 2.5, res = 300, units = "in")
plot.new()
print(p)
mtext("B)",side = 3,line = 2,adj = 0,cex = 1.5,padj = -1)
dev.off()
