source(fn-uvmrplot.R)
# forestplot for all mr method 
# reading results and format data ----------------------------------------------

res = fread(paste0(rdsf_personal,"results/results_bin.csv"))

res_sf = fread(paste0(rdsf_personal,"results/results_sf_bin.csv"))

# draw -------------------------------------------------------------------------
# Options: 
# "Urate (CKDGen)" "eGFR (CKDGen)"  "Urate (UKB)"  "SBP (UKB)"    "DBP (UKB)"

dev.off()
p = uvmr_plot(dat = res,
              exp = "Urate (CKDGen)",
              out = c("eGFR (CKDGen)", "SBP (UKB)", "DBP (UKB)"),
              line_number = 2,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in urate",
              x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.15,0.2),
              intervals = c(-0.2,0.2))

p

pdf(paste0(rdsf_personal,"results/urate on each outcome forestplot.pdf"),width = 16, height = 5)
print(p)
dev.off()

p = uvmr_plot(dat = res,
              exp = "SBP (UKB)",
              out = c("eGFR (CKDGen)", "Urate (CKDGen)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in SBP",
              x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.15,0.2),
              intervals = c(-0.2,0.2))

p

pdf(paste0(rdsf_personal,"results/sbp on each outcome forestplot.pdf"),width = 16, height = 3.5)
print(p)
dev.off()

p = uvmr_plot(dat = res,
              exp = "DBP (UKB)",
              out = c("eGFR (CKDGen)", "Urate (CKDGen)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in DBP",
              x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.15,0.2),
              intervals = c(-0.2,0.2))

p

pdf(paste0(rdsf_personal,"results/dbp on each outcome forestplot.pdf"),width = 16, height = 3.5)
print(p)
dev.off()

p = uvmr_plot(dat = res,
              exp = "eGFR (CKDGen)",
              out = c("SBP (UKB)", "DBP (UKB)", "Urate (CKDGen)"),
              line_number = 2,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in log(eGFR)",
              x_ticks = c(-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1),
              intervals = c(-0.25,0.1))

p

pdf(paste0(rdsf_personal,"results/egfr on each outcome forestplot.pdf"),width = 16, height = 5)
print(p)
dev.off()
     
# additional figures -----------------------------------------------------------

sup_results = fread(paste0(rdsf_personal,"results/sup_results.csv"))

sup_sf_results = fread(paste0(rdsf_personal,"results/sup_sf_results.csv"))

p = uvmr_plot(dat = rbind(sup_results,res %>% filter(exposure %in% c("SBP (UKB)","DBP (UKB)")) %>% filter(outcome == "eGFR (CKDGen)")),
              exp = "SBP (UKB)",
              out = c("eGFR (CKDGen)","Serum creatinine (eGFRcrea)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in log(eGFR)",
              x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2),
              intervals = c(-0.2,0.2))

pdf(paste0(rdsf_personal,"results/compare SBP on eGFR forestplot.pdf"),width = 16, height = 5)
print(p)
dev.off()

p = uvmr_plot(dat = rbind(sup_results,res %>% filter(exposure %in% c("SBP (UKB)","DBP (UKB)")) %>% filter(outcome == "eGFR (CKDGen)")),
              exp = "DBP (UKB)",
              out = c("eGFR (CKDGen)","Serum creatinine (eGFRcrea)"),
              line_number = 1,
              xlabel = "Beta (with 95% CI) for each continuous outcome per SD unit change in log(eGFR)",
              x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2),
              intervals = c(-0.25,0.1))

pdf(paste0(rdsf_personal,"results/compare DBP on eGFR forestplot.pdf"),width = 16, height = 5)
print(p)
dev.off()

