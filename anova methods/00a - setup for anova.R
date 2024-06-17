#' Setup functions and data in a fresh R session

# load functions
source("01 - distributions.R")
source("02 - anova simulations.R")

# read in summary statistics from prior teleneuropsych data
dat <- readRDS("Summary Statistics.RDS")
