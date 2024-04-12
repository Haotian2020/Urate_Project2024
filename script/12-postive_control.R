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

urate_gout_mr <- uvmr("Urate (UKB)","ebi-a-GCST001790",ncase = 2115,ncontrol = 67259)

sbp_stroke_mr <- uvmr("SBP (UKB)","stroke",ncase = 40585,ncontrol = 406111)

dbp_stroke_mr <- uvmr("DBP (UKB)","stroke",ncase = 40585,ncontrol = 406111)

# Save the results -------------------------------------------------------------

pos_res <- rbind(urate_gout_mr[[1]],sbp_stroke_mr[[1]],dbp_stroke_mr[[1]],
                 urate_gout_mr[[4]],sbp_stroke_mr[[4]],dbp_stroke_mr[[4]]) %>% generate_odds_ratios()

pos_pleio <- rbind(urate_gout_mr[[2]],sbp_stroke_mr[[2]],dbp_stroke_mr[[2]],
                   urate_gout_mr[[5]],sbp_stroke_mr[[5]],dbp_stroke_mr[[5]]) %>% generate_odds_ratios()
  
pos_hetero <- rbind(urate_gout_mr[[3]],sbp_stroke_mr[[3]],dbp_stroke_mr[[3]],
                    urate_gout_mr[[6]],sbp_stroke_mr[[6]],dbp_stroke_mr[[6]])%>% generate_odds_ratios()

write.table(pos_res,file = paste0(rdsf_personal,"results/pos_res.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(pos_pleio,file = paste0(rdsf_personal,"results/pos_pleio.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(pos_hetero,file = paste0(rdsf_personal,"results/pos_hetero.csv"),
            sep= ',', row.names = F,col.names= T)

# Plot -------------------------------------------------------------------------

mydata = fread(paste0(rdsf_personal,"results/pos_res.csv")) %>% subset(exposure == "Urate (UKB)"&type == "Ori")%>%drawbinaryformat()
p1 = drawbinaryfigure(mydata,"urate","gout",c(1,4,8,12,16),c(1,16))

mydata = fread(paste0(rdsf_personal,"results/pos_res.csv"))%>%subset(exposure == "SBP (UKB)"&type == "Ori")%>%drawbinaryformat()
p2 = drawbinaryfigure(mydata,"SBP","stroke",c(1,1.5,2,2.5,3),c(1,3))

mydata = fread(paste0(rdsf_personal,"results/pos_res.csv"))%>%subset(exposure == "DBP (UKB)"&type == "Ori")%>%drawbinaryformat()
p3 = drawbinaryfigure(mydata,"DBP","stroke",c(1,1.5,2,2.5,3),c(1,3))

dev.off()
pdf(paste0(rdsf_personal,"results/urate on gout forestplot.pdf"),width = 17, height = 2.5)
plot.new()
mtext("A)",side = 3,line = 2,adj = 0, cex = 1.5,padj = 0)
print(p1)
dev.off()
pdf(paste0(rdsf_personal,"results/sbp on stroke forestplot.pdf"),width = 17, height = 2.5)
plot.new()
print(p2)
mtext("B)",side = 3,line = 2,adj = 0,cex = 1.5,padj = 0)
dev.off()
pdf(paste0(rdsf_personal,"results/dbp on stroke forestplot.pdf"),width = 17, height = 2.5)
plot.new()
print(p3)
mtext("C)",side = 3,line = 2,adj = 0,cex = 1.5,padj = 0)
dev.off()
