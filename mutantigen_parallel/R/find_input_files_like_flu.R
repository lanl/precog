# Determining which input files are in the ballpark of the original MutAntiGen parameters.
## AC Murph
## Date: Apr 2026
library(dplyr)
library(parallel)
library(doParallel)
setwd(here::here())

# Parameter ranges for the current mutantigen lhs:
param_ranges <- list(
  initialNs = c(500000, 42000000),
  demeAmplitudes = c(0, 0.2),
  lambdaAntigenic = c(8.57143e-05, 0.002571429),
  meanAntigenicSize = c(1.2e-3, 1.2e-1),
  lambda = c(0.0095, 4.08),
  mutCost = c(8e-4, 8e-2),
  beta = c(0.1428, 2.25),
  nu = c(1/4, 1/14),
  epsilon_mut = c(0.5, 1.5),
  initialI_prop = c(0.0001, .001)
)
# lhs_df$initialI <- lhs_df$initialNs * lhs_df$initialI_prop
# lhs_df$externalMigration <- lhs_df$initialNs* 0.00001
# lhs_df$epsilon <- (.16/40000000) * lhs_df$initialNs * lhs_df$epsilon_mut

# Parameter values within the original parameters_load.yml:
initialNs = 40000000
demeAmplitudes = 0
lambdaAntigenic = 0.00075
meanAntigenicSize = 0.012
lambda = 0.10
mutCost = 0.008 
beta = 0.5627
nu = 0.25
epsilon = 0.16
initialI = 7406
externalMigration = 200

# Create 10% above and below ranges based on original values
pm10 <- function(x) c(max(0, x * 0.9), x * 1.1)

param_ranges_orig <- list(
  initialNs = pm10(initialNs),
  demeAmplitudes = c(0, 0.1), # I just chose these since otherwise it'd be identically zero.
  lambdaAntigenic = pm10(lambdaAntigenic),
  meanAntigenicSize = pm10(meanAntigenicSize),
  lambda = pm10(lambda),
  mutCost = pm10(mutCost),
  beta = pm10(beta),
  nu = pm10(nu),
  epsilon = pm10(epsilon),
  initialI = pm10(initialI),
  externalMigration = pm10(externalMigration)
)

# With these parameter ranges, let's check all the input files.
# files to inspect
files <- list.files(
  path = "input_files",
  pattern = "^parameters_load_updated_[0-9]+\\.yml$",
  full.names = TRUE
)

# function to extract a scalar numeric value from YAML entry
extract_scalar <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (length(x) == 0) return(NA_real_)
  as.numeric(x[1])
}

# check one file against param_ranges_orig
check_one_file <- function(file, param_ranges) {
  y <- yaml::read_yaml(file)
  
  do.call(rbind, lapply(names(param_ranges), function(param) {
    value <- extract_scalar(y[[param]])
    rng <- param_ranges[[param]]
    
    data.frame(
      file = basename(file),
      parameter = param,
      value = value,
      range_min = rng[1],
      range_max = rng[2],
      in_range = !is.na(value) && value >= rng[1] && value <= rng[2],
      present_in_file = !is.null(y[[param]]),
      stringsAsFactors = FALSE
    )
  }))
}

sockettype <- "PSOCK"
ncores <- 90
cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
file_data <- foreach(file = files,
                  .verbose = TRUE,
                  .errorhandling = "pass",
                  .combine = "rbind"
                ) %dopar% {
                    check_one_file(file, param_ranges_orig)
}
parallel::stopCluster(cl)

file_summary <- file_data %>%
  group_by(file) %>%
  summarise(
    all_params_in_range = all(in_range & present_in_file),
    n_failed = sum(!(in_range & present_in_file))
  )
file_summary


