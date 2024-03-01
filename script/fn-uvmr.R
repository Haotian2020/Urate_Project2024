uvmr <- function(exposure, outcome) {
  
  # Load instruments -----------------------------------------------------------
  
  instruments <- data.table::fread("data/instruments_all.txt", data.table = FALSE)
  
  ## Create exposure dataset ---------------------------------------------------
  
  exp <- TwoSampleMR::format_data(dat = instruments[instruments$id==exposure,],
                                  type = "exposure",
                                  phenotype_col = "exposure")

}

