# function for formatting meta-analysis from two MR results
# e.g. results from split-sample method

format_meta <- function(mr1, mr2) {
  mr_bin = c()
  for(each in mr1$method){
    print(each)
    mr1sub = subset(mr1,method == each)
    mr2sub = subset(mr2,method == each)
    print("do meta-analysis for two MR results")
    meta_mr = metagen(TE = c(mr1sub$b, mr2sub$b), seTE = c(mr1sub$se, mr2sub$se))
    print("writing Meta Results")
    combmr = data.frame(rbind(
      c("META","META",
        paste0(mr1sub$outcome," + ",mr2sub$outcome),
        paste0(mr1sub$exposure," + ",mr2sub$exposure),
        each,
        paste0(mr1sub$nsnp,' & ',mr2sub$nsnp),
        meta_mr$TE.random,
        meta_mr$seTE.random,
        meta_mr$pval.random)))
    
    colnames(combmr) = c("id.exposure","id.outcome",
                         "outcome","exposure",
                         "method",
                         "nsnp","b","se","pval")
    combmr$b = as.numeric(combmr$b)
    combmr$se = as.numeric(combmr$se)
    combmr$pval = as.numeric(combmr$pval)
    mr_bin = rbind(mr_bin,combmr)}
  return(mr_bin)
}


format_mvmrmeta <- function(mr1, mr2,exp1,exp2) {
  mr_bin = c()
  mr1exp1 = mr1$result[stringr::str_detect(string = mr1$result$exposure, pattern = exp1),]
  mr2exp1 = mr2$result[stringr::str_detect(string = mr2$result$exposure, pattern = exp1),]
  
  mr1exp2 = mr1$result[stringr::str_detect(string = mr1$result$exposure, pattern = exp2),]
  mr2exp2 = mr2$result[stringr::str_detect(string = mr2$result$exposure, pattern = exp2),]
  
  print(paste0("do meta-analysis for MVMR results for ",exp1))
  meta_mrexp1 = metagen(TE = c(mr1exp1$b, mr2exp1$b), seTE = c(mr1exp1$se, mr2exp1$se))
  print(paste0("do meta-analysis for MVMR results for ",exp2))
  meta_mrexp2 = metagen(TE = c(mr1exp2$b, mr2exp2$b), seTE = c(mr1exp2$se, mr2exp2$se))
  
  print("writing Meta Results")
  combmr = data.frame(rbind(
    c("META","META", # for exp1
      paste0(mr1exp1$outcome," + ",mr2exp1$outcome),
      paste0(mr1exp1$exposure," + ",mr2exp1$exposure),
      "MVMR",
      paste0(mr1exp1$nsnp,' & ',mr2exp1$nsnp),
      meta_mrexp1$TE.random,
      meta_mrexp1$seTE.random,
      meta_mrexp1$pval.random),
    c("META","META", # for exp2
      paste0(mr1exp2$outcome," + ",mr2exp2$outcome),
      paste0(mr1exp2$exposure," + ",mr2exp2$exposure),
      "MVMR",
      paste0(mr1exp2$nsnp,' & ',mr2exp2$nsnp),
      meta_mrexp2$TE.random,
      meta_mrexp2$seTE.random,
      meta_mrexp2$pval.random)))
  
  colnames(combmr) = c("id.exposure","id.outcome",
                       "outcome","exposure",
                       "method",
                       "nsnp","b","se","pval")
  combmr$b = as.numeric(combmr$b)
  combmr$se = as.numeric(combmr$se)
  combmr$pval = as.numeric(combmr$pval)
  mr_bin = combmr
  return(mr_bin)
}
