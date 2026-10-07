# Shared setup for SIR space-filling experiments
# This script contains all common initialization code used by all four experiment parts
# Author: AC Murph
# Date: Feb 2026

# Load required libraries
library(parallel)
library(doParallel)
library(doSNOW)
library(ggplot2)
library(patchwork)
library(Rcpp)

# Compile and load C++ functions (only once per session)
# sourceCpp() creates the R wrapper functions directly, so we don't need RcppExports.R
cpp_file <- here::here("src", "OSFD.cpp")
if (!exists("approx_gen") || !is.function(approx_gen)) {
  cat("Compiling C++ functions from", cpp_file, "...\n")
  # Add inst/include to the include path for nanoflann.hpp
  Sys.setenv(PKG_CPPFLAGS = paste0("-I", here::here("inst/include")))
  Rcpp::sourceCpp(cpp_file)
  cat("C++ functions loaded successfully\n")
}

# Load helper functions (sourceCpp already loaded C++ functions, so skip RcppExports.R)
source(here::here("R", "source_helpers.R"))
source_helpers_result <- source_helpers()

# Remove RcppExports.R wrappers if they were loaded (they won't work without the package)
# The sourceCpp() versions are the correct ones to use
if (exists("RcppExports.R") && "RcppExports.R" %in% basename(source_helpers_result$files)) {
  cat("Note: RcppExports.R wrapper functions replaced by sourceCpp() versions\n")
}

# Set working directory
setwd(here::here())

# Define all parameters (nsuccesses will be passed as argument to each script)
# nsuccesses should be set by the calling script before sourcing this
# cand_batch should be set in scripts that use OSFD (parts 3 & 4)
p = 3
q = 2
mc.cores = 99
CAND_size = 5000000
pkg_location = here::here()
alpha_bounds = c(0.001, 100)
reproduction_number_bounds = c(1.001, 100)
pit_bounds = c(1, 50)
piv_bounds = c(0.001, 1)
s0_bounds = c(0.95, 0.999)

# Feasibility checking function is already loaded via source_helpers()

cat("\n=== SIR Experiment Setup Complete ===\n")
cat(sprintf("mc.cores: %d\n", mc.cores))
cat(sprintf("Parameters: p=%d, q=%d\n", p, q))
cat("=====================================\n\n")
