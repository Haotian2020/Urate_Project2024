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
  return(dat)
}