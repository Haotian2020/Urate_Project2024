# forestplot for all mr method 
# reading results and format data ----------------------------------------------

res = fread(paste0(urate_github,"result/results_bin.csv"))

# change the exposure and outcome names for figures ----------------------------

replacement_rules <- c("urate_clean" = "Urate UKB", "sbp_clean" = "SBP UKB", "dbp_clean" = "DBP UKB")

for (pattern in names(replacement_rules)) {
  res$exposure <- sub(pattern, replacement_rules[pattern], res$exposure)
  res$outcome <- sub(pattern, replacement_rules[pattern], res$outcome)
}

# draw -------------------------------------------------------------------------
# Options: 
# "Urate CKDGen" "eGFR CKDGen"  "Urate UKB"  "SBP UKB"    "DBP UKB"

p = uvmr_plot(res,"Urate CKDGen",c("eGFR CKDGen", "SBP UKB", "DBP UKB"))

dev.off()
png(paste0(github,"figures/"),width = 15000, height = 5000,res = 1000)
print(p)
dev.off()
