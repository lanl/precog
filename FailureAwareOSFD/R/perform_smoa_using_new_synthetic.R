# Create the real-data embeddings used to do the experiments in the supplement
## Author: LJ Beesley + AC Murph
## Date: January 2025
library(ggplot2)
library(data.table)
library(plyr)
library(gridExtra)
library(lubridate)
library(parallel)
library(doParallel)
library(doSNOW)
library(grid)
library(plotly)
library(GGally)
library(dplyr)
library(tidyr)
library(mgcv)
library(collapse)
theme_set(theme_bw())

# Load helper functions
source(here::here("R", "source_helpers.R"))
source_helpers()

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(task_id)){
  task_id = 1
}

setwd(here::here())

data_path = here::here("data", "raw_data")
output_path = here::here("data", "embeddings_gam_real")

####################
### Load work manifest for flattened parallelization
####################
manifest_file <- here::here("data", "work_manifest.csv")
if (!file.exists(manifest_file)) {
  stop("Work manifest not found. Run R/create_work_manifest.R first.")
}

manifest <- read.csv(manifest_file)
cat(sprintf("Loaded manifest with %d total work items\n", nrow(manifest)))

# Number of SLURM array tasks (nodes)
n_nodes <- 6

# Divide manifest into equal chunks for each node
total_items <- nrow(manifest)
items_per_node <- ceiling(total_items / n_nodes)

# Calculate which rows this task should process
start_idx <- (task_id - 1) * items_per_node + 1
end_idx <- min(task_id * items_per_node, total_items)

my_work <- manifest[start_idx:end_idx, ]

cat(sprintf("\n=== SLURM Task %d ===\n", task_id))
cat(sprintf("Processing work items %d to %d (%d items)\n", start_idx, end_idx, nrow(my_work)))
cat(sprintf("Diseases in this chunk: %s\n", paste(unique(my_work$disease), collapse=", ")))

####################
### Prepare sMOA ###
####################
k <- 5
h <- 4
closest <- 2000
ncores <- 99
sockettype <- "PSOCK"

# Load pre-computed synthetic embeddings
embeddings_file <- here::here('data', 'synthetic_embeddings.RData')
if (!file.exists(embeddings_file)) {
  stop("Synthetic embeddings not found. Run R/prepare_synthetic_embeddings.R first.")
}
load(embeddings_file)
print(paste("Loaded embeddings from", embeddings_file))
print(paste("ISFD dimensions: X =", paste(dim(embed_mat_X_ISFD), collapse=" x "),
            ", y =", paste(dim(embed_mat_y_ISFD), collapse=" x ")))
print(paste("OSFD dimensions: X =", paste(dim(embed_mat_X_OSFD), collapse=" x "),
            ", y =", paste(dim(embed_mat_y_OSFD), collapse=" x ")))

####################
### Pre-load all disease data for this chunk (avoid redundant I/O)
####################
unique_diseases <- unique(my_work$disease)
cat(sprintf("\nPre-loading %d unique diseases for this chunk...\n", length(unique_diseases)))
disease_data <- list()
for (disease in unique_diseases) {
  cat(sprintf("  Loading %s...\n", disease))
  disease_data[[disease]] <- readRDS(file.path(data_path, paste0(disease, ".RDS")))
  cat(sprintf("    Loaded %d time series\n", length(disease_data[[disease]])))
}
cat("All disease data loaded into memory\n")

## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoSNOW(cl)

# Create progress bar
pb <- txtProgressBar(max = nrow(my_work), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

cat(sprintf("Processing %d time series using %d cores\n", nrow(my_work), ncores))

# Export all necessary data to parallel workers
train_data <- foreach(i = 1:nrow(my_work),
      .errorhandling = "pass",
      .verbose = T,
      .options.snow = opts,
      .export = c('disease_data', 'my_work', 'embed_mat_X_OSFD', 'embed_mat_y_OSFD',
                  'embed_mat_X_ISFD', 'embed_mat_y_ISFD', 'data_path', 'output_path',
                  'k', 'h', 'closest'),
      .packages = c('dplyr', 'tsfeatures', 'timeDate', 'lubridate','mgcv','zoo','forecast','collapse'))%dopar%{

        # Get disease and time series index for this work item
        eval_key <- my_work$disease[i]
        j <- my_work$ts_index[i]

        # Get pre-loaded disease data (no redundant disk I/O)
        list_of_lists <- disease_data[[eval_key]]

        item <- list_of_lists[[j]]
        ts_test = list_of_lists[[j]]
        fcst_indices <- 7:(length(ts_test$ts)-h)

        # Rolling window question?   What regime should we associate with ts_test

        lastobs = NULL
        full_ret_mat = NULL
        if(!file.exists(paste0(output_path,"/embed_",eval_key, '_',j,"_full_ret_mat.csv"))){
          if((max(ts_test$ts) > 10 & ts_test$ts_scale == 'counts') | (max(ts_test$ts)>1e-3 & ts_test$ts_scale == 'proportion')){
            for(l in fcst_indices){
              ## get info_packet and trim
              info_packet <- ts_test
              info_packet$ts <- info_packet$ts[1:l]
              data_till_now = data.frame(value = info_packet$ts, t = 1:l)
              if(nrow(data_till_now)>50){
                data_till_now = data_till_now[(nrow(data_till_now)-50):nrow(data_till_now),]
              }

              # Check if last k and future h values are all zero (variation <= 1e-8)
              last_k_vals <- tail(info_packet$ts, k)
              future_h_vals <- ts_test$ts[(l+1):(l+h)]
              if(var(last_k_vals) <= 1e-8 && var(future_h_vals) <= 1e-8) {
                next  # Skip this forecast
              }

              #### light smoothing and differencing and get last k
              data_till_now_smoothed      <- gam(value~ s(t,k=pmax(round(nrow(data_till_now)/2))),data=data_till_now)$fitted.values

              #### Classify regime for entire smoothed series
              regime_classifications <- classify_entire_timeseries(data_till_now_smoothed, snippet_length = k)

              #### could use gam smoother
              to_match_in_moa             <- tail((data_till_now_smoothed),k)
              mu_to_match = mean(to_match_in_moa)
              sd_to_match = (sd(to_match_in_moa)+1e-8)
              to_match_in_moa = (to_match_in_moa - mu_to_match) / sd_to_match

              # Get the regime for the current position (last regime in the classification)
              current_regime <- tail(regime_classifications$regime, 1) # This carries the last value forward to the end of the timeseries, so tail should be okay here.  It is the classification of the snipped starting k spaces ago.

              y_temp = ts_test$ts[(l+1):(l+h)]

              #### make MOA w new synthetic
              ts = to_match_in_moa
              ret_mat <- data.frame(h = 1:h)
              dist_to_test = rowSums(abs(embed_mat_X_OSFD %r-% tail(ts,k))) 
              min_dist = order(dist_to_test)[1:closest]
              # min_dist <- sort(dist_to_test,index.return = TRUE)$ix[1:closest]
              # x1 = embed_mat_X[min_dist[1],]

              ret_mat$moa_OSFD <- pmax(0,apply(embed_mat_y_OSFD[min_dist,],2,median))*sd_to_match + mu_to_match
              ret_mat$min_dist_OSFD = min(dist_to_test) #min_dist
              # summary(dist_to_test)
              rm('dist_to_test')

              #### make MOA w old synthetic
              ts = to_match_in_moa
              dist_to_test = rowSums(abs(embed_mat_X_ISFD %r-% tail(ts,k))) 
              min_dist = order(dist_to_test)[1:closest]
              # min_dist <- sort(dist_to_test,index.return = TRUE)$ix[1:closest]
              # x2 = X[min_dist[1],]

              ret_mat$moa_ISFD <- pmax(0,apply(embed_mat_y_ISFD[min_dist,],2,median))*sd_to_match + mu_to_match
              ret_mat$min_dist_ISFD = min(dist_to_test) #min_dist
              # summary(dist_to_test)
              rm('dist_to_test')

              # Revert scaling prior to recording obs:
              # to_match_in_moa             <- tail((data_till_now_smoothed),k)
              ret_mat$obs = data_till_now$value[length(data_till_now$value)]
              ret_mat$truth = y_temp
              ret_mat$target_end_date = info_packet$ts_dates[l]
              ret_mat$regime_of_real_data = current_regime

              full_ret_mat = rbind(full_ret_mat, ret_mat)
            }
            write.csv(full_ret_mat, file = paste0(output_path,"/embed_",eval_key, '_',j,"_full_ret_mat.csv"), quote = F, row.names = F)
          }

        }

        # Return summary info
        list(
          disease = eval_key,
          ts_index = j,
          processed = !is.null(full_ret_mat)
        )
}
close(pb)
stopCluster(cl)

# Print summary
cat("\n=== Processing Complete ===\n")
successful <- sum(sapply(train_data, function(x) !inherits(x, "try-error")))
cat(sprintf("Successfully processed: %d / %d time series\n", successful, nrow(my_work)))
cat(sprintf("Task %d complete\n", task_id))
