# function to plot mr results

uvmr_plot <- function(dat, exp, out, line_number, xlabel) {
  
  mydata  <- data.frame(dat) %>%
     subset(.,exposure == exp, outcome == out)
  
# remove id.exposure and id.outcome columns ------------------------------------
  
  if(!is.na(mydata$id.exposure & !is.na(mydata$id.outcome))){
    mydata = select(mydata, -c("id.exposure","id.outcome"))%>%generate_odds_ratios()
  }
  
# format data ------------------------------------------------------------------
  
  mydata$beta <-
    format(round(mydata$b, digits = 3),
           nsmall = 3)
  mydata$CI <-
    paste0(format(round(mydata$lo_ci, digits = 3),
                  nsmall = 3),
           ", ",
           format(round(mydata$up_ci, digits = 3),
                  nsmall = 3))
  mydata$pvalue <-
    format(round(mydata$pval, digits = 3))
  mydata$method <-
    factor(
      mydata$method,
      levels = c(
        "Inverse variance weighted","Weighted median","Steiger Filtering",
        "Weighted mode","MR Egger", "Simple mode"))
  
  sorted_index <- order(mydata$method)
  mydata = mydata[sorted_index,]
  
# make the table shown with the figure -----------------------------------------
  
  for (i in c(seq(2, nrow(mydata), 5),
              seq(3, nrow(mydata), 5),
              seq(4, nrow(mydata), 5),
              seq(5, nrow(mydata), 5))) {
    mydata$outcome[i] <- NA
    mydata$exposure[i] <- NA
    mydata$nsnp[i] <- NA
  }
  
  mydata = data.frame(mydata)
  
  tabletext <- cbind(
    c("Expsoure", as.character(mydata[, 'exposure'])),
    c("Outcome", as.character(mydata[, 'outcome'])),
    c("Approach", as.character(mydata[, 'method'])),
    c('Number of SNPs', as.character(mydata[, 'nsnp'])),
    c("Beta", as.character(mydata[, 'beta'])),
    c("95% CI", as.character(mydata[, 'CI'])),
    c("p-value", as.character(mydata[, 'pvalue'])))
  
  dev.off()
  p <- forestplot(
    tabletext,
    graph.pos = 4,
    mean = as.numeric(rbind(NA, cbind(mydata[, 'beta']))),
    lower = as.numeric(rbind(NA, cbind(mydata[, 'lo_ci']))),
    upper = as.numeric(rbind(NA, cbind(mydata[, 'up_ci']))),
    new_page = F,
    psignif = 0.05,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    hrzl_lines = list(),
    
    boxsize = 0.15,
    line.margin = 0.1,
    lty.ci = 1,
    col = fpColors(box = "black", lines = "darkgray"),
    lwd.ci = 1,
    ci.vertices = T,
    ci.vertices.height = 0.15,
    graphwidth = unit(150, "mm"),
    is.summary = c(T, rep(F, nrow(tabletext))),
    colgap = unit (5, "mm"),
    zero = 0.0,
    xticks = c(-0.2,-0.1,0,0.1,0.2,0.3),
    clip = c(-0.2,0.3),
    xlab = paste0(""))
  
   return(p)
}
