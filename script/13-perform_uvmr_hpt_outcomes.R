# Script for performing MR of urate on early, late -onset and overall hypertension

source("fn-uvmr.R")
source("fn-binaryplot.R")

# read instruments files -------------------------------------------------------

all_int <- fread(paste0(rdsf_personal,"/data/format_data/all_instruments.csv")) %>% data.frame()

# Create empty results datasets ------------------------------------------------

hpt_res <- NULL 
hpt_hetero <- NULL
hpt_pleio <- NULL

# Perform MR -------------------------------------------------------------------

outcomes = c("early50","late60","hpt")

for(i in outcomes){
  if(i == "early50"){
    x = 6934
    y = 329146
    urate_hpt_mr_bin <- uvmr("exurate_sd",outcome = i,ncase = x,ncontrol = y)
  }else if(i == "late60"){
    x = 95583
    y = 329146
    urate_hpt_mr_bin <- uvmr("exurate_sd",outcome = i,ncase = x,ncontrol = y)
  }else if(i == "hpt"){
    x = 133680
    y = 329146
    urate_hpt_mr_bin <- uvmr("exurate_sd",outcome = i,ncase = x,ncontrol = y)
  }
  hpt_res = rbind(hpt_res,urate_hpt_mr_bin[[1]],urate_hpt_mr_bin[[4]])
  hpt_pleio = rbind(hpt_pleio,urate_hpt_mr_bin[[2]],urate_hpt_mr_bin[[5]])
  hpt_hetero = rbind(hpt_hetero,urate_hpt_mr_bin[[3]],urate_hpt_mr_bin[[6]])
}

# Save the results -------------------------------------------------------------

write.table(hpt_res,file = paste0(rdsf_personal,"results/hpt_res.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(hpt_pleio,file = paste0(rdsf_personal,"results/hpt_pleio.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(hpt_hetero,file = paste0(rdsf_personal,"results/hpt_hetero.csv"),
            sep= ',', row.names = F,col.names= T)

# Plot -------------------------------------------------------------------------

mydata = fread(paste0(rdsf_personal,"results/hpt_res.csv")) %>% subset(type == "Ori")%>% generate_odds_ratios()
mydata = drawbinaryformat(mydata)

for (i in c(seq(2, nrow(mydata), 5),
            seq(3, nrow(mydata), 5),
            seq(4, nrow(mydata), 5),
            seq(5, nrow(mydata), 5))) {
  mydata$outcome[i] <- NA
  mydata$exposure[i] <- NA
  mydata$nsnp[i] <- NA}
mydata$exposure[-1] <- NA
mydata = data.frame(mydata)

tabletext <- cbind(
  c("Expsoure", as.character(mydata[, 'exposure'])),
  c("Outcome", as.character(mydata[, 'outcome'])),
  c("Approach", as.character(mydata[, 'method'])),
  c('Number of SNPs', as.character(mydata[, 'nsnp'])),
  c("OR", as.character(mydata[, 'OR'])),
  c("95% CI", as.character(mydata[, 'CI'])),
  c("p-value", as.character(mydata[, 'pvalue'])))

p <- forestplot(
  tabletext,
  graph.pos = 4,
  mean = as.numeric(rbind(NA, cbind(mydata[, 'or']))),
  lower = as.numeric(rbind(NA, cbind(mydata[, 'or_lci95']))),
  upper = as.numeric(rbind(NA, cbind(mydata[, 'or_uci95']))),
  new_page = F,
  psignif = 0.05,
  xlog = F,
  txt_gp = fpTxtGp(
    label = gpar(cex = 1),
    ticks = gpar(cex = 1),
    xlab = gpar(cex = 1),
    title = gpar(cex = 1)
  ),
  hrzl_lines = list("7" = gpar(lty = 1, lwd = 1, col = "black"),
                    "12" = gpar(lty = 1, lwd = 1, col = "black")
  ),
  boxsize = 0.3,
  line.margin = 0.2,
  lty.ci = 1,
  col = fpColors(box = "black", lines = "darkgray"),
  lwd.ci = 2,
  ci.vertices = T,
  ci.vertices.height = 0.15,
  graphwidth = unit(150, "mm"),
  is.summary = c(T, rep(F, nrow(tabletext))),
  colgap = unit (5, "mm"),
  zero = 1.0,
  xticks = c(0.8,1,1.2,1.4,1.6),
  clip = c(0.8,1.6),
  xlab = paste0("OR (with 95% CI) for all types of hypertension per SD unit change in urate"))
p

dev.off()
pdf(paste0(rdsf_personal,"results/exurate on all hpt forestplot.pdf"),width = 17, height = 5)
print(p)
dev.off()