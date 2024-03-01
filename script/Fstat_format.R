F_statistic = function(exposure_data){
  exposure_data$R = get_r_from_bsen(exposure_data$beta.exposure,exposure_data$se.exposure,exposure_data$samplesize.exposure)
  exposure_data$Rsq = exposure_data$R^2
  exposure_data$F_stat = (exposure_data$samplesize.exposure-2)*((exposure_data$Rsq)/(1-exposure_data$Rsq))
  exposure_data
  
  return(exposure_data)
}

# read data from excel
all_exp = read_excel(sheet = 3, path = paste0(github,"results/Supplementary Tables_Urate_final.xlsx"),skip =1)
all_exp$R = get_r_from_bsen(all_exp$beta.exposure,all_exp$se.exposure,all_exp$samplesize.exposure)
all_exp$Rsq = all_exp$R^2
all_exp$F_stat = (all_exp$samplesize.exposure-2)*((all_exp$Rsq)/(1-all_exp$Rsq))
all_exp

write.table(all_exp,file = paste0(github,"results/all_exp.csv"),
            sep= ',', row.names = F,col.names= T)
