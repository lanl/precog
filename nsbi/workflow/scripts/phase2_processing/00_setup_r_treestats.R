#!/usr/bin/env Rscript

#' One-time setup: Install treestats and its remaining CRAN dependencies
#' 
#' This script runs ONCE before any parallel batch jobs to avoid race conditions.
#' It installs the 5 packages not available in conda-forge:
#'   - pso
#'   - subplex
#'   - treebalance  
#'   - DDD
#'   - treestats
#'
#' All other dependencies (31 packages) are pre-installed from conda.

cat("========================================\n")
cat("Setting up R treestats environment\n")
cat("========================================\n\n")

# Packages that must be installed from CRAN (not available in conda)
cran_packages <- c("pso", "subplex", "treebalance", "DDD", "treestats")

cat("Checking for packages not available in conda...\n")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from CRAN...\n", pkg))
    install.packages(
      pkg, 
      repos = "https://cloud.r-project.org",
      Ncpus = 4,
      quiet = FALSE
    )
    
    # Verify installation
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Failed to install %s", pkg))
    }
    cat(sprintf("✓ %s installed successfully\n\n", pkg))
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

cat("\n========================================\n")
cat("Testing treestats package...\n")
cat("========================================\n")

# Load and test treestats
library(treestats)
library(ape)

# Create a simple test tree
test_tree <- rtree(10)
cat("Created test tree with", Ntip(test_tree), "tips\n")

# Test a simple statistic
colless <- colless(test_tree)
cat("Colless index:", colless, "\n")

cat("\n✓ treestats setup complete!\n")
cat("All packages ready for parallel batch processing.\n")
