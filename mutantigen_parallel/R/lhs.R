##################################################################### 
#' Script for generating the parameters via lhs and for writing
#' the yaml codes to be read in by MutAntiGen
#' Step 1: Set parameter bounds
#' Step 2: Run the LHS and recast to get the parameter alues
#' Step 3: Update any dependent parameters 
#' Step 4: Save
#' Step 5: Set up paths and directories for yaml files
#' Step 6: Loop through the parameter settings

library(lhs)   # For running the parameter lhs
library(yaml)  # For working with yaml files 
library(readr) # For reading data
library(pracma)
library(caret)
library(gbm)
set.seed(0225)
setwd(paste0(this.path::here(), "/../"))

#### Step 1: Determine the parameter settings #### 

# Define the number of samples
n_samples <- 4000

# Define the parameter ranges
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

#### Step 2:  Generate the Latin Hypercube Sample ####
lhs_samples <- randomLHS(n_samples, length(param_ranges))

# Scale the LHS samples to the defined parameter ranges

# For (initialNs) - uniform across the log scale 
range       <- param_ranges[[1]]
log_samples <- (log10(range[2]) - log10(range[1])) * lhs_samples[,1] + log10(range[1])
log_samples <- 10^log_samples

# For rest of the samples
scaled_samples <- t(apply(lhs_samples, 1, function(x) {
  sapply(2:length(param_ranges), function(i) {
    range <- param_ranges[[i]]
    range[1] + x[i] * (range[2] - range[1])
  })
}))

scaled_samples <- cbind(log_samples, scaled_samples)


# Convert the scaled samples into a data frame for readability
colnames(scaled_samples) <- names(param_ranges)
lhs_df <- as.data.frame(scaled_samples)


#### Step 3: Parameters that have to be adjusted based on others ####
lhs_df$thresholdAntigenicSize <- lhs_df$meanAntigenicSize
lhs_df$initialI <- lhs_df$initialNs * lhs_df$initialI_prop
lhs_df$externalMigration <- lhs_df$initialNs* 0.00001
lhs_df$epsilon <- (.16/40000000) * lhs_df$initialNs * lhs_df$epsilon_mut
# lhs_df <- lhs_df[,-which(colnames(lhs_df) == "initialI_prop")]

#lhs_model = readRDS(file = "gbm_f1_step1.RDS")
# lhs_model = readRDS(file = "gbm_fbeta_step2.RDS")

#completion_predictions = which(predict(object=lhs_model, newdata=as.matrix(lhs_df))=="X1")
# murph note: I'm taking this out for this run.
modified_lhs = lhs_df#[completion_predictions,]
nrow(modified_lhs)

# Preview the sampled data
# modified_lhs %>% 
#   pivot_longer(names_to = "param", values_to = "values", cols = everything()) %>% 
#   ggplot(aes(values)) + geom_histogram() + facet_wrap(~param, scales = "free") + 
#   scale_x_log10()

#### Step 4: Save to a CSV if needed ####
write.csv(modified_lhs, "lhs_samples.csv", row.names = FALSE)

#### Step 5: Now update the yaml files ########################################

# Function to enforce correct data types for the YAML content
update_parameter <- function(param_name, new_value, yaml_content) {
  if (param_name %in% c("burnin", "endDay", "printStep", "tipSamplingStartDay", "tipSamplingEndDay",
                        "tipSamplesPerDeme", "demeCount", "initialI", "initialDeme",
                        "antigenicEvoStartDay", "hostImmuneHistorySampleCount", "fitSampleCount",
                        "printFitSamplesStep")) {
    yaml_content[[param_name]] <- as.integer(new_value)
  } else if (param_name %in% c("tipSamplingRate", "treeProportion", "yearsFromMK",
                               "birthRate", "deathRate", "beta", "nu", "betweenDemePro",
                               "externalMigration", "immunityLoss", "initialPrT", "backgroundDistance",
                               "muPhenotype","smithConversion", "homologousImmunity", "initialTraitA", "meanStep",
                               "sdStep", "lambda", "mutCost", "probLethal", "epsilon", "epsilonSlope",
                               "lambdaAntigenic", "meanAntigenicSize", "antigenicGammaShape",
                               "thresholdAntigenicSize", "cleanUpDistance", "demoNoiseScaler")) {
    yaml_content[[param_name]] <- as.numeric(new_value)
  } else if (param_name %in% c("tipSamplingProportional", "repeatSim", "immunityReconstruction",
                               "memoryProfiling", "pcaSamples", "detailedOutput",
                               "restartFromCheckpoint", "swapDemography", "transcendental",
                               "backgroundImmunity", "mut2D")) {
    yaml_content[[param_name]] <- as.logical(new_value)
  } else if (param_name %in% c("demeNames")) {
    # Ensure `demeNames` is a list of strings
    yaml_content[[param_name]] <- as.list(strsplit(new_value, ",\\s*")[[1]])
  } else if (param_name %in% c("initialNs")) { 
    # Ensure arrays are correctly typed (integer or double)
    if (!is.character(new_value)) {
      new_value <- as.character(new_value)
    }
    yaml_content[[param_name]] <- as.list(as.integer(strsplit(new_value, ",\\s*")[[1]]))
  } else if (param_name %in% c("demeAmplitudes")) {
    yaml_content[[param_name]] <- as.list(as.numeric(new_value)) 
  }
  else if (param_name %in% c("phenotypeSpace")) {
    # `phenotypeSpace` is expected as a plain string
    yaml_content[[param_name]] <- as.character(new_value)
  } else {
    # For other parameters, update directly
    yaml_content[[param_name]] <- new_value
  }
  return(yaml_content)
}

# Set up paths and directories
yaml_file_path <- "parameters_load.yml"
csv_file_path <- "lhs_samples.csv"

# Read the YAML file into a list
yaml_content <- yaml.load_file(yaml_file_path)

# Read the CSV file containing the new parameters
new_parameters <- read_csv(csv_file_path)
print(paste0("The modified LHS sample gave ", nrow(new_parameters), 
             " parameter locations that are predicted to complete." ))

# Create a directory to save the modified YAML files
output_dir <- "input_files/"

#### Step 6: Loop to update the YAML content using new parameters
for (i in 1:nrow(new_parameters)) {
  # Create a copy of the original YAML content for this set of parameters
  updated_yaml <- yaml_content
  
  # Update the YAML content using the current row of the CSV
  for (param_name in colnames(new_parameters)) {
    new_value <- new_parameters[[param_name]][i]
    updated_yaml <- update_parameter(param_name, new_value, updated_yaml)
  }
  
  # Force formatting of content that doesn't change so not in parameter list - for now
  value <- yaml_content[["demeBaselines"]]
  updated_yaml[["demeBaselines"]] <- as.list(as.numeric(value)) 
  
  value <- yaml_content[["demeOffsets"]]
  updated_yaml[["demeOffsets"]] <- as.list(as.numeric(value)) 
  
  value <- yaml_content[["demeNames"]]
  updated_yaml[["demeNames"]] <- as.list(strsplit(value, ",\\s*")[[1]])
  
  # Save the updated YAML file
  output_file <- file.path(output_dir, paste0("parameters_load_updated_", i, ".yml"))
  write_yaml(updated_yaml, output_file)
  print(paste("Finished writing yaml file number", i))
}

cat("Updated YAML files have been saved to:", output_dir, "\n")
