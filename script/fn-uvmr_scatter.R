# function for plotting scatter plot for UVMR results

save_scatter_plot = function(dat, res, sf = 0){
  
  p1 <- mr_scatter_plot(res, dat)
  
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)
  
  res_loo <- mr_leaveoneout(dat)
  p3 <- mr_leaveoneout_plot(res_loo)
  
  res_single <- mr_singlesnp(dat)
  p4 <- mr_funnel_plot(res_single)
  
  if(sf == 0){
  ggsave(p1[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_scatterplot.pdf"))
  ggsave(p2[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_forestplot.pdf"))
  ggsave(p3[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_leaveoneout.pdf"))
  ggsave(p4[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"_funnelplot.pdf"))}
  
  if(sf == 1){
    ggsave(p1[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"sf_scatterplot.pdf"))
    ggsave(p2[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"sf_forestplot.pdf"))
    ggsave(p3[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"sf_leaveoneout.pdf"))
    ggsave(p4[[1]], file = paste0(rdsf_personal,"results/",unique(res$exposure),"_on_",unique(res$outcome),"sf_funnelplot.pdf"))
  }
  
}
