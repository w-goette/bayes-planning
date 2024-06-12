#' Setup functions and data in a fresh R session

# load functions
source("01 - distributions.R")
source("02a - t statistics.R")
source("02b - t posteriors.R")
source("03 - t simulations.R")
source("04 - plot t-test.R")

# read in summary statistics from prior teleneuropsych data
dat <- readRDS("Summary Statistics.RDS")
