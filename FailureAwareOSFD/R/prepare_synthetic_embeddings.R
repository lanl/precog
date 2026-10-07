# Prepare synthetic embeddings from stochastic SIR data
## Author: LJ Beesley + AC Murph
## Date: January 2025
library(parallel)
library(doParallel)
source(here::here("R", "source_helpers.R"))
source_helpers()

####################
### Prepare sMOA ###
####################
k <- 5
h <- 4
ncores <- 90
include_original_smoa_lib = FALSE
sockettype <- "PSOCK"

# Gather all the synthetic data and create embedding matrices
embed_mat_X_OSFD = NULL
embed_mat_y_OSFD = NULL
regime_OSFD = NULL

embed_mat_X_ISFD = NULL
embed_mat_y_ISFD = NULL
regime_ISFD = NULL

if(include_original_smoa_lib){
  cat("Loading original sMOA library data...\n")
  load(here::here('data', 'synthetic_simidx_1_num_curves_18387_orig.RData'))
  cat("Creating embedding matrix from original data...\n")
  embed_mat                       <- create_embed_matrix(sim_ts,h=(h),k=(k))
  X                               <- embed_mat[[1]]
  y                               <- embed_mat[[2]]
  cat("Original data loaded: X has", nrow(X), "rows,", ncol(X), "columns\n")
}


# Get all CSV files from data/stochastic_sir/
cat("\n=== Starting CSV file processing ===\n")
sir_data_dir <- here::here('data', 'stochastic_sir')
if (dir.exists(sir_data_dir)) {
  cat("Directory found:", sir_data_dir, "\n")
  all_files <- list.files(sir_data_dir, pattern = "\\.csv$", full.names = TRUE)
  cat("Total CSV files found:", length(all_files), "\n")

  # Separate ISFD and OSFD files
  isfd_files <- all_files[grepl("_ISFD\\.csv$", all_files)]
  osfd_files <- all_files[grepl("_OSFD\\.csv$", all_files)]
  cat("ISFD files:", length(isfd_files), "\n")
  cat("OSFD files:", length(osfd_files), "\n")

  # Set up parallel cluster for file reading
  cat("\nSetting up parallel cluster with", ncores, "cores...\n")
  cl_read <- parallel::makeCluster(spec = ncores, type = sockettype)
  registerDoParallel(cl_read)
  cat("Cluster ready.\n")

  # Process ISFD files in parallel
  if (length(isfd_files) > 0) {
    cat("\n=== Processing ISFD files ===\n")
    cat("Starting parallel processing of", length(isfd_files), "ISFD files...\n")
    flush.console()

    # Progress tracking
    progress_file <- tempfile()

    isfd_results <- foreach(file_name = isfd_files, .combine = 'c', .multicombine = TRUE,
                            .options.snow = list(progress = function(n) {
                              if (n %% 100 == 0) {
                                cat(sprintf("  Processed %d/%d ISFD files (%.1f%%)\n",
                                           n, length(isfd_files), 100*n/length(isfd_files)))
                                flush.console()
                              }
                            })) %dopar% {
      data <- read.csv(file_name)
      time_cols <- grep("^t[0-9]+$", colnames(data))
      if (length(time_cols) >= (k + h)) {
        # Extract as matrices to ensure consistent structure
        X_data <- as.matrix(data[, time_cols[1:k]])
        y_data <- as.matrix(data[, time_cols[(k+1):(k+h)]])
        list(list(X = X_data,
                  y = y_data,
                  regime = data[,1]))
      } else {
        list(NULL)
      }
    }
    cat("Parallel processing complete. Combining results...\n")
    flush.console()

    # Combine results
    isfd_results <- isfd_results[!sapply(isfd_results, is.null)]
    cat("Valid ISFD results:", length(isfd_results), "\n")
    if (length(isfd_results) > 0) {
      cat("Combining ISFD X matrices...\n")
      embed_mat_X_ISFD <- do.call(rbind, lapply(isfd_results, function(x) x$X))
      cat("Combining ISFD y matrices...\n")
      embed_mat_y_ISFD <- do.call(rbind, lapply(isfd_results, function(x) x$y))
      cat("Combining ISFD regime labels...\n")
      regime_ISFD <- do.call(c, lapply(isfd_results, function(x) x$regime))
      cat("ISFD combined dimensions: X =", paste(dim(embed_mat_X_ISFD), collapse=" x "), "\n")

      # Include the original sMOA library (potentially)
      if(include_original_smoa_lib){
        cat("Merging original sMOA library with ISFD data...\n")
        # Convert original data to matrix with matching column names
        X_mat <- as.matrix(X)
        colnames(X_mat) <- colnames(embed_mat_X_ISFD)
        embed_mat_X_ISFD = rbind(embed_mat_X_ISFD, X_mat)

        y_mat <- as.matrix(y)
        colnames(y_mat) <- colnames(embed_mat_y_ISFD)
        embed_mat_y_ISFD = rbind(embed_mat_y_ISFD, y_mat)

        regime_ISFD = c(regime_ISFD, rep(NA, times = nrow(X)))
        cat("After merging: ISFD X =", paste(dim(embed_mat_X_ISFD), collapse=" x "), "\n")
      }

      # Scale X rows by their row-mean and row-std
      cat("Scaling ISFD X by row statistics...\n")
      X_row_means <- rowMeans(embed_mat_X_ISFD)
      X_row_sds <- apply(embed_mat_X_ISFD, 1, sd)
      embed_mat_X_ISFD <- (embed_mat_X_ISFD - X_row_means) / (X_row_sds+1e-6)

      # Scale y using X row statistics
      cat("Scaling ISFD y using X statistics...\n")
      embed_mat_y_ISFD <- (embed_mat_y_ISFD - X_row_means) / (X_row_sds+1e-6)
      cat("ISFD scaling complete.\n")
    }
  }

  # Process OSFD files in parallel
  if (length(osfd_files) > 0) {
    cat("\n=== Processing OSFD files ===\n")
    cat("Starting parallel processing of", length(osfd_files), "OSFD files...\n")
    flush.console()

    osfd_results <- foreach(file_name = osfd_files, .combine = 'c', .multicombine = TRUE,
                            .options.snow = list(progress = function(n) {
                              if (n %% 100 == 0) {
                                cat(sprintf("  Processed %d/%d OSFD files (%.1f%%)\n",
                                           n, length(osfd_files), 100*n/length(osfd_files)))
                                flush.console()
                              }
                            })) %dopar% {
      data <- read.csv(file_name)
      time_cols <- grep("^t[0-9]+$", colnames(data))
      if (length(time_cols) >= (k + h)) {
        # Extract as matrices to ensure consistent structure
        X_data <- as.matrix(data[, time_cols[1:k]])
        y_data <- as.matrix(data[, time_cols[(k+1):(k+h)]])
        list(list(X = X_data,
                  y = y_data,
                  regime = data[,1]))
      } else {
        list(NULL)
      }
    }
    cat("Parallel processing complete. Combining results...\n")
    flush.console()

    # Combine results
    osfd_results <- osfd_results[!sapply(osfd_results, is.null)]
    cat("Valid OSFD results:", length(osfd_results), "\n")
    if (length(osfd_results) > 0) {
      cat("Combining OSFD X matrices...\n")
      embed_mat_X_OSFD <- do.call(rbind, lapply(osfd_results, function(x) x$X))
      cat("Combining OSFD y matrices...\n")
      embed_mat_y_OSFD <- do.call(rbind, lapply(osfd_results, function(x) x$y))
      cat("Combining OSFD regime labels...\n")
      regime_OSFD <- do.call(c, lapply(osfd_results, function(x) x$regime))
      cat("OSFD combined dimensions: X =", paste(dim(embed_mat_X_OSFD), collapse=" x "), "\n")

      if(include_original_smoa_lib){
        cat("Merging original sMOA library with OSFD data...\n")
        # Convert original data to matrix with matching column names
        X_mat <- as.matrix(X)
        colnames(X_mat) <- colnames(embed_mat_X_OSFD)
        embed_mat_X_OSFD = rbind(embed_mat_X_OSFD, X_mat)

        y_mat <- as.matrix(y)
        colnames(y_mat) <- colnames(embed_mat_y_OSFD)
        embed_mat_y_OSFD = rbind(embed_mat_y_OSFD, y_mat)

        regime_OSFD = c(regime_OSFD, rep(NA, times = nrow(X)))
        cat("After merging: OSFD X =", paste(dim(embed_mat_X_OSFD), collapse=" x "), "\n")
      }

      # Scale X rows by their row-mean and row-std
      cat("Scaling OSFD X by row statistics...\n")
      X_row_means <- rowMeans(embed_mat_X_OSFD)
      X_row_sds <- apply(embed_mat_X_OSFD, 1, sd)
      embed_mat_X_OSFD <- (embed_mat_X_OSFD - X_row_means) / (X_row_sds+1e-6)

      # Scale y using X row statistics
      cat("Scaling OSFD y using X statistics...\n")
      embed_mat_y_OSFD <- (embed_mat_y_OSFD - X_row_means) / (X_row_sds+1e-6)
      cat("OSFD scaling complete.\n")
    }
  }

  # Stop the cluster
  cat("\nStopping parallel cluster...\n")
  parallel::stopCluster(cl_read)
  cat("Cluster stopped.\n")
}

cat("\n=== Final Results ===\n")
cat("ISFD Breakdown:\n")
print(table(regime_ISFD))
cat("\nOSFD Breakdown:\n")
print(table(regime_OSFD))

# Save the embeddings to data directory
cat("\nSaving results...\n")
output_dir <- here::here('data')
output_file <- file.path(output_dir, 'synthetic_embeddings.RData')
save(embed_mat_X_ISFD, embed_mat_y_ISFD, regime_ISFD,
     embed_mat_X_OSFD, embed_mat_y_OSFD, regime_OSFD,
     file = output_file)

cat("\n=== SUCCESS ===\n")
cat("Saved embeddings to:", output_file, "\n")
cat("ISFD dimensions: X =", paste(dim(embed_mat_X_ISFD), collapse=" x "),
            ", y =", paste(dim(embed_mat_y_ISFD), collapse=" x "), "\n")
cat("OSFD dimensions: X =", paste(dim(embed_mat_X_OSFD), collapse=" x "),
            ", y =", paste(dim(embed_mat_y_OSFD), collapse=" x "), "\n")
cat("Script complete!\n")
