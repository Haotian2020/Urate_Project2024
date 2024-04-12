# function for plotting scatter plot for UVMR results

save_scatter_plot = function(dat,res){
  p1 <- mr_scatter_plot(res, dat)
  
  ggsave(p1[[1]], file = paste0(rdsf_personal,"results/",unique(mr$exposure),"_on_",unique(mr$outcome),"_scatterplot.pdf"), width = 7, height = 7)
  
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)
  
  res_loo <- mr_leaveoneout(dat)
  p3 <- mr_leaveoneout_plot(res_loo)
  
  res_single <- mr_singlesnp(dat)
  p4 <- mr_funnel_plot(res_single)
  
}
