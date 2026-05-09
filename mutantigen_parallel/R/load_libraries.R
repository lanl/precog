## Dave Osthus
## 6-4-25
## 4 quadrants paper: libraries


# ################################################################################
# # load_packages.R
# # Automatically installs and loads required R packages
# ################################################################################
# 
# # List of required packages
# packages <- c(
#   "ggplot2",
#   "data.table",
#   "mgcv",
#   "lightgbm",
#   "reshape2",
#   "gridExtra",
#   "grid",
#   "lubridate",
#   "zoo",
#   "Matrix",
#   "pls",
# # "forecast",
#   "foreach",
#   "parallel",
#   "doParallel",
#   "RANN",
#   "FNN",
#   "RColorBrewer",
#   "colorspace",
#   "scales",
#   "dplyr",
#   "ggnewscale",
#   "tools",
#   "cowplot",
#   "uwot",
#   "fastcluster",
#   "ggrepel"
# )
# 
# # Function to check if a package is installed; if not, install it, then load it
# install_if_missing <- function(pkg) {
#   # Check if package is available without loading it
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)  # Install if missing
#   }
#   library(pkg, character.only = TRUE)  # Load the package using its name as a string
# }
# 
# # Apply the function to each package in the list
# # `invisible()` suppresses printing the result of lapply to the console
# invisible(lapply(packages, install_if_missing))





# List of required packages
packages <- c(
  "ggplot2",
  "data.table",
  "mgcv",
  "lightgbm",
  "reshape2",
  "gridExtra",
  "grid",
  "lubridate",
  "zoo",
  "Matrix",
  "pls",
  "foreach",
  "parallel",
  "doParallel",
  "RANN",
  "FNN",
  "RColorBrewer",
  "colorspace",
  "dplyr",
  "tools",
  "cowplot",
  "uwot",
  "ggrepel",
  "scales",
  "doSNOW",
  "utils",
  "MASS",
  "fields",
  "pathviewr",
  "lhs",
  "ks",
  "hexbin",
  "shadowtext",
  "truncnorm",
  "tibble",
  "tidyr"
)

# Function to install and load a package with messages
install_and_load <- function(pkg) {
  #cat("Processing package:", pkg, "\n")
  
  # Try to install if not available
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  -> Package not found. Attempting to install...\n")
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = T)
      #cat("  -> Installation successful.\n")
    }, error = function(e) {
      cat("  -> Installation failed:", conditionMessage(e), "\n")
      return(FALSE)
    })
  }
  
  # Try to load the package
  success <- tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) {
    cat("  -> Loading failed:", conditionMessage(e), "\n")
    FALSE
  })
  
  if (success) {
    #cat("  -> Loaded successfully.\n\n")
  } else {
    cat("  -> Could not load package.\n\n")
  }
}

# Run through each package one by one
for (pkg in packages) {
  install_and_load(pkg)
}
