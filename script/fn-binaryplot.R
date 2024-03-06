drawbinaryformat= function(dat){
  
  dat = dat %>% select(-c("id.exposure","id.outcome")) %>% split_outcome() %>%  generate_odds_ratios()
  
  dat$OR <-
    format(round(dat$or, digits = 3),
           nsmall = 3)
  dat$CI <-
    paste0(format(round(dat$or_lci95, digits = 3),
                  nsmall = 3),
           ", ",
           format(round(dat$or_uci95, digits = 3),
                  nsmall = 3))
  dat$pvalue <-
    format(round(dat$pval, digits = 3))
  dat$method <-
    factor(
      dat$method,
      levels = c(
        "Inverse variance weighted","Weighted median",
        "Weighted mode","MR Egger", "Simple mode"))
  
  dat$outcome <-
    factor(
      dat$outcome,
      levels = c("Gout","Stroke",
                 "Hypertension (UKB)",
                 "Early-onset hypertension (UKB)",
                 "Late-onset hypertension (UKB)"
        ))
  
  sorted_index <- order(dat$outcome,dat$method)
  
  dat = dat[sorted_index,]
  
  for (i in c(seq(2, nrow(dat), 5),
              seq(3, nrow(dat), 5),
              seq(4, nrow(dat), 5),
              seq(5, nrow(dat), 5))) {
    dat$outcome[i] <- NA
    dat$exposure[i] <- NA
    dat$nsnp[i] <- NA}
  dat = data.frame(dat)
  return(dat)
  
  }

drawbinaryfigure = function(dat,exposure,outcome,xseq,xminxmax){
  tabletext <- cbind(
    c("Expsoure", as.character(dat[, 'exposure'])),
    c("Outcome", as.character(dat[, 'outcome'])),
    c("Approach", as.character(dat[, 'method'])),
    c('Number of SNPs', as.character(dat[, 'nsnp'])),
    c("OR", as.character(dat[, 'OR'])),
    c("95% CI", as.character(dat[, 'CI'])),
    c("p-value", as.character(dat[, 'pvalue'])))
  
  p <- forestplot(
    tabletext,
    graph.pos = 4,
    mean = as.numeric(rbind(NA, cbind(dat[, 'or']))),
    lower = as.numeric(rbind(NA, cbind(dat[, 'or_lci95']))),
    upper = as.numeric(rbind(NA, cbind(dat[, 'or_uci95']))),
    new_page = F,
    psignif = 0.05,
    xlog = F,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    hrzl_lines = list(
      
    ),
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
    zero = 1.0,
    xticks = xseq,
    clip = xminxmax,
    xlab = paste0("OR (with 95% CI) for ",outcome," per SD unit change in ",exposure))
  return(p)
}
