# function for lines in the figures

generate_lines <- function(line_number) {
  lines <- list()
  
  for (i in seq_len(line_number)) {
    key <- as.character(i * 5 + 2)
    lines[[key]] <- gpar(lty = 1, lwd = 1, col = "black")
  }
  
  return(lines)
}

# function to plot UVMR results

uvmr_plot <- function(dat, exp, out, line_number, xlabel, x_ticks, intervals){
  # remove id.exposure and id.outcome columns ------------------------------------
  
  mydata  <- data.frame(dat) %>%
    dplyr::select(-c("id.exposure", "id.outcome")) %>%
    subset(exposure %in% exp & outcome %in% out) %>%
    generate_odds_ratios()
  
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
        "Inverse variance weighted",
        "MR Egger",
        "Weighted median",
        "Weighted mode",
        "Simple mode",
        "Steiger Filtering"
      )
    )
  
  mydata$exposure <-
    factor(
      mydata$exposure,
      levels = c(
        "Urate (CKDGen)",
        "Urate (UKB)",
        "SBP (UKB)",
        "DBP (UKB)",
        "PP (UKB)",
        "BMI",
        "eGFR (CKDGen2019)",
        "eGFR (CKDGen2016)"
      )
    )

  mydata$outcome <-
    factor(
      mydata$outcome,
      levels = c(
        "Urate (CKDGen)",
        "Urate (UKB)",
        "SBP (UKB)",
        "DBP (UKB)",
        "PP (UKB)",
        "BMI",
        "eGFR (CKDGen2019)",
        "eGFR (CKDGen2016)"
      )
    )
  
  mydata <- mydata %>%
    mutate(method = if_else(method == "Inverse variance weighted", "IVW", method))
 
  if ("type" %in% names(mydata) && all(c("Ori", "Steiger") %in% mydata$type)) {
    mydata$type <- factor(
      mydata$type,
      levels = c("Ori", "Steiger")
    )
  }
  
  sorted_index <- order(mydata$outcome,mydata$type,mydata$method)
  mydata = mydata[sorted_index, ]
  
  # make the table shown with the figure -----------------------------------------
  
  for (i in c(seq(2, nrow(mydata), 5),
              seq(3, nrow(mydata), 5),
              seq(4, nrow(mydata), 5),
              seq(5, nrow(mydata), 5))) {
    mydata$outcome[i] <- NA
    mydata$exposure[i] <- NA
    mydata$nsnp[i] <- NA
  }
  mydata[-1, "exposure"] <- NA
  mydata = data.frame(mydata)
  
  tabletext <- cbind(
    c("Expsoure", as.character(mydata[, 'exposure'])),
    c("Outcome", as.character(mydata[, 'outcome'])),
    c("Approach", as.character(mydata[, 'method'])),
    c('Number of SNPs', as.character(mydata[, 'nsnp'])),
    c("Beta", as.character(mydata[, 'beta'])),
    c("95% CI", as.character(mydata[, 'CI'])),
    c("p-value", as.character(mydata[, 'pvalue']))
  )
  
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
    hrzl_lines = generate_lines(line_number),
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
    zero = 0.0,
    xticks = x_ticks,
    clip = intervals,
    xlab = paste0(xlabel)
  )
  
  return(p)
}
