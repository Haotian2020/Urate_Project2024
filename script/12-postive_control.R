# Script for doing positive control MR analyses for the outcome GWAS conducted in UKB
# urate on gout
# bp on stroke

source("fn-binaryplot.R")

# read instruments files -------------------------------------------------------

all_int <- read_tsv(paste0(rdsf_personal,"/data/format_data/all_instruments.csv")) %>% data.frame()

# Create empty results datasets ------------------------------------------------

pos_res <- NULL
pos_hetero <- NULL
pos_pleio <- NULL

pos_res_sf <- NULL
pos_hetero_sf <- NULL
pos_pleio_sf <- NULL

# Perform MR -------------------------------------------------------------------
mr_results <- list(
  uvmr("urate", "ebi-a-GCST001790", ncase = 2115, ncontrol = 67259),
  uvmr("sbp", "stroke", ncase = 40585, ncontrol = 406111),
  uvmr("dbp", "stroke", ncase = 40585, ncontrol = 406111),
  uvmr("pp", "stroke", ncase = 40585, ncontrol = 406111),
  uvmr("egfr_sd", "ckdmvp", ncase = 63705, ncontrol = 368937)
)

combine_results <- function(idx1, idx2) {
  do.call(rbind, lapply(mr_results, function(res) rbind(res[[idx1]], res[[idx2]])))
}

# merge results
pos_res <- generate_odds_ratios(combine_results(1, 4))
pos_pleio <- combine_results(2, 5)
pos_hetero <- combine_results(3, 6)

# save files
write.table(pos_res, file = paste0(rdsf_personal, "results/pos_res.csv"),
            sep = ',', row.names = FALSE, col.names = TRUE)
write.table(pos_pleio, file = paste0(rdsf_personal, "results/pos_pleio.csv"),
            sep = ',', row.names = FALSE, col.names = TRUE)
write.table(pos_hetero, file = paste0(rdsf_personal, "results/pos_hetero.csv"),
            sep = ',', row.names = FALSE, col.names = TRUE)

# Plot -------------------------------------------------------------------------

all_data <- fread(paste0(rdsf_personal, "results/pos_res.csv"))

# define parameters in the plot
plot_list <- list(
  list(exposure = "Urate (UKB)", lab = "urate", outcome = "gout", xvals = c(1,4,8,12), ylim = c(1,12), letter = "A", filename = "urate on gout forestplot.pdf"),
  list(exposure = "SBP (UKB)", lab = "SBP", outcome = "stroke", xvals = c(1,1.5,2,2.5), ylim = c(1,2.5), letter = "B", filename = "sbp on stroke forestplot.pdf"),
  list(exposure = "DBP (UKB)", lab = "DBP", outcome = "stroke", xvals = c(1,1.5,2,2.5), ylim = c(1,2.5), letter = "C", filename = "dbp on stroke forestplot.pdf"),
  list(exposure = "PP (UKB)", lab = "PP", outcome = "stroke", xvals = c(1,1.5,2,2.5), ylim = c(1,2.5), letter = "D", filename = "pp on stroke forestplot.pdf"),
  list(exposure = "eGFR", lab = "PP", outcome = "stroke", xvals = c(1,1.5,2,2.5), ylim = c(1,2.5), letter = "D", filename = "pp on stroke forestplot.pdf"),
)


for (pinfo in plot_list) {
  mydata <- subset(all_data, exposure == pinfo$exposure & type == "Ori")
  p <- drawbinaryfigure(mydata, pinfo$lab, pinfo$outcome, pinfo$xvals, pinfo$ylim)
  pdf(paste0(rdsf_personal, "results/", pinfo$filename), width = 17, height = 2.5)
  plot.new()
  mtext(paste0(pinfo$letter, ")"), side = 3, line = 2, adj = 0, cex = 1.5, padj = 0)
  print(p)
  dev.off()
}
