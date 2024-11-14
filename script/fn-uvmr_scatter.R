# function for plotting scatter plot for UVMR results

save_scatter_plot = function(dat,res){
  
  p1 <- mr_scatter_plot(res, dat)
  ggsave(p1[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_scatterplot.pdf"))
  
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)
  ggsave(p2[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_forestplot.pdf"))
  
  res_loo <- mr_leaveoneout(dat)
  p3 <- mr_leaveoneout_plot(res_loo)
  ggsave(p3[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_leaveoneout.pdf"))
  
  res_single <- mr_singlesnp(dat)
  p4 <- mr_funnel_plot(res_single)
  ggsave(p4[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_funnelplot.pdf"))

}
