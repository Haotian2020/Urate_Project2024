# Script for plotting UVMR figures

source(fn-uvmrplot.R)

# reading results and format data ----------------------------------------------
res = fread(paste0(rdsf_personal,"results/results_bin.csv"))

res_sf = fread(paste0(rdsf_personal,"results/results_sf_bin.csv"))

# draw -------------------------------------------------------------------------
# Options for exposure: 
# "Urate (CKDGen)" "eGFR (CKDGen2019)"  "Urate (UKB)"  "SBP (UKB)"   "DBP (UKB)" "PP (UKB)"

dev.off()
p = uvmr_plot(dat = res,
              exp = "Urate (CKDGen)",
              out = c("eGFR (CKDGen2019)", "SBP (UKB)", "DBP (UKB)", "PP (UKB)"),
              line_number = 3,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in urate",
              x_ticks = c(-0.2,-0.1,0,0.1,0.2),
              intervals = c(-0.2,0.2))

pdf(paste0(rdsf_personal,"results/exurate on each outcome forestplot.pdf"),width = 18, height = 6)
plot.new()
mtext("A)",side = 3,line = 2,adj = 0, cex = 1.5,padj = 0)
print(p)
dev.off()

p = uvmr_plot(dat = res,
              exp = "SBP (UKB)",
              out = c("eGFR (CKDGen2019)", "Urate (CKDGen)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in SBP",
              x_ticks = c(-0.2,-0.1,0,0.1,0.2),
              intervals = c(-0.2,0.2))

p

pdf(paste0(rdsf_personal,"results/sbp on each outcome forestplot.pdf"),width = 18, height = 3.5)
plot.new()
mtext("B)",side = 3,line = 2,adj = 0, cex = 1.5,padj = 0)
print(p)
dev.off()

p = uvmr_plot(dat = res,
              exp = "DBP (UKB)",
              out = c("eGFR (CKDGen2019)", "Urate (CKDGen)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in DBP",
              x_ticks = c(-0.2,-0.1,0,0.1,0.2),
              intervals = c(-0.2,0.2))

p

pdf(paste0(rdsf_personal,"results/dbp on each outcome forestplot.pdf"),width = 18, height = 3.5)
plot.new()
mtext("C)",side = 3,line = 2,adj = 0, cex = 1.5,padj = 0)
print(p)
dev.off()

# add pp plot
p = uvmr_plot(dat = res,
              exp = "PP (UKB)",
              out = c("eGFR (CKDGen2019)", "Urate (CKDGen)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in PP",
              x_ticks = c(-0.2,-0.1,0,0.1,0.2),
              intervals = c(-0.2,0.2))

p

pdf(paste0(rdsf_personal,"results/pp on each outcome forestplot.pdf"),width = 18, height = 3.5)
plot.new()
mtext("D)",side = 3,line = 2,adj = 0, cex = 1.5,padj = 0)
print(p)
dev.off()

p = uvmr_plot(dat = res,
              exp = "eGFR (CKDGen2019)",
              out = c("SBP (UKB)", "DBP (UKB)", "PP (UKB)", "Urate (CKDGen)"),
              line_number = 3,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in log(eGFR)",
              x_ticks = c(-0.2,-0.1,0,0.1,0.2),
              intervals = c(-0.2,0.1))

p

pdf(paste0(rdsf_personal,"results/egfr on each outcome forestplot.pdf"),width = 18, height = 6)
plot.new()
mtext("E)",side = 3,line = 2,adj = 0, cex = 1.5,padj = 0)
print(p)
dev.off()