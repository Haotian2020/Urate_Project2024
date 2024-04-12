# function for cleaning continuous trait, removing outliers located outside of 4SD

singleconti_clean <- function(dat,var){
  print(paste0("total data set No.row is ",nrow(dat)))
  dat = select(dat,c("eid",all_of(var))) %>% na.omit()
  print(paste0("after removing NA, No.row is ",nrow(dat)))
  print(colnames(dat))
  dat = merge(dat, link, by = "eid")
  print(paste0("after merging with link, No.row is ",nrow(dat)))
  dat = dat[,c(3,3,2)]
  print(head(dat))
  colnames(dat) = c("FID","IID",var)
  var_mean = mean(dat[[var]])
  var_sd = sd(dat[[var]])
  # print(paste0("mean of ",var," is ",var_mean))
  # print(paste0("sd of ",var," is ",var_sd ))
  dat_clean = subset(dat,dat[[3]]<=(var_mean+4*var_sd) & dat[[3]]>=(var_mean - 4*var_sd))
  print(paste0("the number of outliers were removed ", nrow(dat) - nrow(dat_clean)))
  print(paste0("mean of ",var," is ",mean(dat_clean[[3]])))
  print(paste0("sd of ",var," is ",sd(dat_clean[[3]])))
  print(paste0("cleaned data No.row is ",nrow(dat_clean)))
  return(dat_clean)
}