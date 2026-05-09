# Run the MutAntiGen Simulation in Parallel using R for file organization and custom terminal commands
## Author: AC Murph
## Date: Dec 2024

library(tidyverse)
library(parallel)
library(doParallel)
library(this.path)

# Reset the time logfiles.
# Specify the directory
dir_path <- "time_logfiles"

setwd(paste0(this.path::here(), '/../'))

# List all files in the directory
files <- list.files(dir_path, full.names = TRUE)

# Delete all files
file.remove(files)

# java -cp dist/MutAntiGen.jar:/lib/colt.jar Mutantigen "parameters_load.yml" 1 
ncores <- 99
sockettype <- "PSOCK"

# I'm keeping this for this version so the amount of data is around the same.
number_of_reps = 20
good_input_files = as_tibble(read.csv('antipeak.csv'))

cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:(nrow(good_input_files)*number_of_reps),
                  .verbose = T)%dopar%{
		    start_time  <- Sys.time()
		    
		    # num_of_yaml_row = (i-1) %% nrow(good_input_files) + 1
		    # num_of_yaml = good_input_files$runnum[num_of_yaml_row]

        name_of_yaml = paste0("parameters_load_orig.yml")
		    max_runtime <- 10 * 60 * 60  # 10 hours in seconds
		    full_command <- paste0("timeout ", max_runtime, "s java -cp dist/MutAntiGen.jar:/lib/colt.jar Mutantigen ", name_of_yaml, " ", i)
        # full_command = paste0("java -cp dist/MutAntiGen.jar:/lib/colt.jar Mutantigen ", name_of_yaml, " ", i)
                    
				system(full_command)
		    
		    log_file = paste0('time_logfiles/job_', i, '.txt')
		    # Log the start time
		    cat("Script started at:", start_time, "\n", file = log_file)

 		    # End time
		    end_time <- Sys.time()

		    # Calculate the time difference
		    time_diff <- difftime(end_time, start_time, units = "mins")

		    # Log the end time and duration
		    cat("Script ended at:", end_time, "\n", file = log_file, append = TRUE)
 		    cat("Time taken (mins):", time_diff, "\n", file = log_file, append = TRUE)		    
                  }
stopCluster(cl)
