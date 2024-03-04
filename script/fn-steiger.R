# Function to use Steiger filtering
# Source: https://github.com/MRCIEU/TwoSampleMR/blob/master/R/steiger_filtering.R

steiger_function <- function(dat) {
  print("Please confirm that units of both exposure and outcome are in SD")
  dat$units.exposure <- "SD"
  dat$units.outcome <- "SD"
  #run SF analysis
  dat <- steiger_filtering(dat)
  print(summary(dat$steiger_dir))
  #call subset of data passed steiger filtering
  mr_sf = mr(subset(dat, steiger_dir))
  mr_hetero_sf <- mr_heterogeneity(subset(dat, steiger_dir))
  mr_pleio_sf <- mr_pleiotropy_test(subset(dat, steiger_dir))
  
  # store all the results
  results <- NULL
  results <- list(mr_sf, mr_hetero_sf, mr_pleio_sf)
  
  return(results)
}