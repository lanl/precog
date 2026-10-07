# depends: extract_mutantigen_10_outputs
#' Run MutAntiGen Simulation and Extract Output Summary
#'
#' @description
#' Wrapper function to run a single MutAntiGen simulation with specified parameters,
#' wait for completion, and extract 11-dimensional output summary.
#'
#' @param x Numeric vector of length 10 with the following parameters (in order):
#'   1. initialNs - Initial population size (500000 to 42000000)
#'   2. demeAmplitudes - Amplitude of seasonal forcing (0 to 0.2)
#'   3. lambdaAntigenic - Rate of antigenic evolution (8.57143e-05 to 0.002571429)
#'   4. meanAntigenicSize - Mean size of antigenic jumps (1.2e-3 to 1.2e-1)
#'   5. lambda - Mutation rate (0.0095 to 4.08)
#'   6. mutCost - Fitness cost of mutations (8e-4 to 8e-2)
#'   7. beta - Transmission rate (0.1428 to 2.25)
#'   8. nu - Recovery rate (calculated as beta/R0, where R0 is sampled from 1.01 to 5.0)
#'   9. epsilon_mut - Mutation effect multiplier (0.5 to 1.5)
#'   10. initialI_prop - Initial infected proportion (0.0001 to 0.001)
#'
#' Note: The input vector x should contain nu directly (already calculated from beta/R0)
#'
#' @param NUM Integer simulation ID for file naming
#' @param mutantigen_dir Path to mutantigen_parallel directory (default: auto-detect)
#' @param max_runtime Maximum runtime in seconds (default: 36000 = 10 hours)
#' @param sim_length_days Simulation length in days (default: NULL uses base YAML value of 3650)
#' @param method Character string identifying the sampling method (e.g., "osfd", "lhs", "random")
#'        Used for timing log filenames to compare different approaches
#' @param timing_log_dir Directory to save timing logs (default: "timing_logs" in mutantigen_dir)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return Numeric vector of length 10 with output summaries, or vector of NAs on failure
#'
#' @details
#' This function:
#' 1. Creates a temporary YAML file with the specified parameters
#' 2. Runs the MutAntiGen Java simulation with timeout
#' 3. Waits for completion
#' 4. Extracts 10 output summaries using extract_mutantigen_10_outputs()
#' 5. Returns NA vector if any step fails
#'
#' The 10 outputs are:
#'  1. Total antigenic types >= 1% of population
#'  2. Unique dominant types across peaks
#'  3. Number of peaks (> 5% of population)
#'  4. Average meanR (entire timeseries)
#'  5. Average meanR (FWHM regions only)
#'  6. Dominant type changes between peaks
#'  7. Max peak height
#'  8. Average unique antigens per peak
#'  9. Total infection burden
#' 10. Range of peak heights (max - min)
#'
#' The function is designed to be used with constrained_osfd_ei_with_retry() for
#' output-space filling design experiments.
#'
#' @export
run_mutantigen_and_extract <- function(
  x,
  NUM,
  mutantigen_dir = NULL,
  max_runtime = 36000,  # 10 hours in seconds
  sim_length_days = NULL,
  method = "unknown",
  timing_log_dir = NULL,
  verbose = TRUE
) {
  # ---------------------------
  # Input validation
  # ---------------------------
  if (length(x) != 10) {
    if (verbose) cat("ERROR: Input vector must have length 10\n")
    return(rep(NA_real_, 10))
  }

  # Auto-detect mutantigen directory if not provided
  if (is.null(mutantigen_dir)) {
    # We're in R/, go up one level to find mutantigen_parallel
    pkg_root <- here::here()
    mutantigen_dir <- file.path(pkg_root, "mutantigen_parallel")
    if (!dir.exists(mutantigen_dir)) {
      if (verbose) cat("ERROR: Cannot find mutantigen_parallel directory at:", mutantigen_dir, "\n")
      return(rep(NA_real_, 11))
    }
  }

  # Store original directory
  original_wd <- getwd()

  # Initialize timing log
  timing <- list()
  timing$NUM <- NUM
  timing$method <- method
  timing$start_time <- Sys.time()
  timing$start_timestamp <- format(timing$start_time, "%Y-%m-%d %H:%M:%S")
  timing$success <- FALSE
  timing$exit_code <- NA
  timing$error_message <- ""

  # Wrap everything in tryCatch to ensure we return NA on any error
  result <- tryCatch({
    # ---------------------------
    # Step 1: Parse input parameters
    # ---------------------------
    param_names <- c("initialNs", "demeAmplitudes", "lambdaAntigenic", "meanAntigenicSize",
                     "lambda", "mutCost", "beta", "nu", "epsilon_mut", "initialI_prop")

    params <- setNames(as.list(x), param_names)

    # Calculate dependent parameters
    params$thresholdAntigenicSize <- params$meanAntigenicSize
    params$initialI <- as.integer(params$initialNs * params$initialI_prop)
    params$externalMigration <- params$initialNs * 0.00001
    params$epsilon <- (0.16 / 40000000) * params$initialNs * params$epsilon_mut
    params$initialPrR <- 0.4  # 40% prior immunity - balance between epidemic dynamics and computational speed

    # Reduce tree sampling moderately to speed up simulations
    params$tipSamplingRate <- 0.5          # Sample every 2 days instead of daily
    params$tipSamplesPerDeme <- 50000      # Cap at 50K samples instead of 100K
    params$treeProportion <- 0.5           # Use 50% of tips for tree reconstruction

    if (verbose) {
      cat("===================================\n")
      cat("Running MutAntiGen simulation NUM =", NUM, "\n")
      cat("Method:", method, "\n")
      cat("===================================\n")
    }

    # Store input parameters for logging
    timing$input_params <- as.list(x)
    names(timing$input_params) <- param_names

    # ---------------------------
    # Step 2: Load base YAML template
    # ---------------------------
    t_yaml_start <- Sys.time()
    setwd(mutantigen_dir)

    base_yaml_path <- "parameters_load.yml"
    if (!file.exists(base_yaml_path)) {
      if (verbose) cat("ERROR: Cannot find base YAML file:", base_yaml_path, "\n")
      return(rep(NA_real_, 10))
    }

    yaml_content <- yaml::yaml.load_file(base_yaml_path)

    # ---------------------------
    # Step 3: Update YAML with parameters
    # ---------------------------
    # Helper function to set correct data types
    update_parameter <- function(param_name, new_value, yaml_content) {
      if (param_name %in% c("initialI", "tipSamplesPerDeme", "diversitySamplingCount")) {
        yaml_content[[param_name]] <- as.integer(new_value)
      } else if (param_name %in% c("beta", "nu", "externalMigration", "lambda", "mutCost",
                                   "epsilon", "lambdaAntigenic", "meanAntigenicSize",
                                   "thresholdAntigenicSize", "tipSamplingRate", "treeProportion", "initialPrR")) {
        yaml_content[[param_name]] <- as.numeric(new_value)
      } else if (param_name %in% c("initialNs")) {
        yaml_content[[param_name]] <- list(as.integer(new_value))
      } else if (param_name %in% c("demeAmplitudes")) {
        yaml_content[[param_name]] <- list(as.numeric(new_value))
      } else {
        yaml_content[[param_name]] <- new_value
      }
      return(yaml_content)
    }

    # Update all parameters
    for (param_name in names(params)) {
      yaml_content <- update_parameter(param_name, params[[param_name]], yaml_content)
    }

    # Override simulation length if specified
    if (!is.null(sim_length_days)) {
      yaml_content[["endDay"]] <- as.integer(sim_length_days)
      # Also adjust tipSamplingEndDay to be slightly before endDay
      yaml_content[["tipSamplingEndDay"]] <- as.integer(sim_length_days - 10)
      if (verbose) cat("Overriding simulation length to", sim_length_days, "days\n")
    }

    # Force formatting of arrays that don't change
    if (!is.null(yaml_content[["demeBaselines"]])) {
      yaml_content[["demeBaselines"]] <- as.list(as.numeric(yaml_content[["demeBaselines"]]))
    }
    if (!is.null(yaml_content[["demeOffsets"]])) {
      yaml_content[["demeOffsets"]] <- as.list(as.numeric(yaml_content[["demeOffsets"]]))
    }
    if (!is.null(yaml_content[["demeNames"]])) {
      if (is.character(yaml_content[["demeNames"]])) {
        yaml_content[["demeNames"]] <- list(yaml_content[["demeNames"]])
      }
    }

    # ---------------------------
    # Step 4: Write temporary YAML file
    # ---------------------------
    yaml_dir <- "input_files"
    if (!dir.exists(yaml_dir)) {
      dir.create(yaml_dir, recursive = TRUE)
    }

    yaml_file <- file.path(yaml_dir, paste0("parameters_load_temp_", NUM, ".yml"))
    yaml::write_yaml(yaml_content, yaml_file)

    timing$yaml_creation_time_sec <- as.numeric(difftime(Sys.time(), t_yaml_start, units = "secs"))

    if (verbose) cat("Created YAML file:", yaml_file, "\n")

    # ---------------------------
    # Step 5: Run MutAntiGen simulation
    # ---------------------------
    # Create log directory for Java output
    log_dir <- "logfiles"
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE)
    }

    java_log <- file.path(log_dir, paste0("java_", NUM, ".log"))

    java_cmd <- paste0(
      "timeout ", max_runtime, "s ",
      "java -Xmx32g -Xms16g -cp dist/MutAntiGen.jar:/lib/colt.jar Mutantigen ",
      yaml_file, " ", NUM,
      " > ", java_log, " 2>&1"  # Redirect stdout and stderr to log file
    )

    if (verbose) {
      cat("Running command:\n  ", java_cmd, "\n")
      cat("Max runtime:", max_runtime / 3600, "hours\n")
      cat("Java output will be saved to:", java_log, "\n")
    }

    t_sim_start <- Sys.time()
    exit_code <- system(java_cmd)
    timing$simulation_time_sec <- as.numeric(difftime(Sys.time(), t_sim_start, units = "secs"))
    timing$exit_code <- exit_code

    if (verbose) {
      cat("Simulation finished in", round(timing$simulation_time_sec, 1), "seconds\n")
      cat("Exit code:", exit_code, "\n")
    }

    # Check for timeout or error
    if (exit_code == 124) {
      if (verbose) cat("ERROR: Simulation timed out after", max_runtime, "seconds\n")
      return(rep(NA_real_, 10))
    }
    if (exit_code != 0) {
      if (verbose) cat("ERROR: Simulation failed with exit code", exit_code, "\n")
      return(rep(NA_real_, 10))
    }

    # ---------------------------
    # Step 6: Check output files exist
    # ---------------------------
    output_dir <- "outfiles"
    required_files <- c(
      paste0("out_", NUM, ".timeseries"),
      paste0("out_", NUM, ".viralFitnessSeries"),
      paste0("out_", NUM, ".tips"),
      paste0("out_", NUM, ".out.simmapAntigenic"),
      paste0("out_", NUM, ".simmapLoad")
    )

    missing_files <- character(0)
    for (f in required_files) {
      if (!file.exists(file.path(output_dir, f))) {
        missing_files <- c(missing_files, f)
      }
    }

    if (length(missing_files) > 0) {
      if (verbose) {
        cat("ERROR: Missing output files:\n")
        cat(paste("  -", missing_files, collapse = "\n"), "\n")
      }
      return(rep(NA_real_, 10))
    }

    if (verbose) cat("All output files found\n")

    # ---------------------------
    # Step 7: Extract output summaries
    # ---------------------------
    if (verbose) cat("Extracting output summaries...\n")

    # Note: extraction function is already loaded via source_helpers()
    t_extract_start <- Sys.time()

    outputs <- extract_mutantigen_10_outputs(
      NUM = NUM,
      directory = output_dir,
      first_month_duration = 0.077,  # 4 weeks in years (28/365)
      peak_threshold_pct = 0.01,
      type_threshold_pct = 0.01
    )

    timing$extraction_time_sec <- as.numeric(difftime(Sys.time(), t_extract_start, units = "secs"))

    # Check for failed extraction (all NA)
    if (all(is.na(outputs))) {
      if (verbose) cat("ERROR: Output extraction returned all NAs\n")
      timing$error_message <- "Extraction returned all NAs"
      return(rep(NA_real_, 10))
    }

    # Mark as successful
    timing$success <- TRUE
    timing$output_values <- as.list(outputs)

    if (verbose) {
      cat("Successfully extracted outputs:\n")
      for (i in seq_along(outputs)) {
        cat(sprintf("  Output %2d: %.4f\n", i, outputs[i]))
      }
      cat("===================================\n")
    }

    return(outputs)

  }, error = function(e) {
    if (verbose) {
      cat("ERROR: Exception occurred:\n")
      cat("  ", conditionMessage(e), "\n")
    }
    timing$error_message <- conditionMessage(e)
    return(rep(NA_real_, 10))
  }, finally = {
    # ---------------------------
    # Finalize timing and save log
    # ---------------------------
    timing$end_time <- Sys.time()
    timing$total_time_sec <- as.numeric(difftime(timing$end_time, timing$start_time, units = "secs"))
    timing$total_time_min <- timing$total_time_sec / 60
    timing$end_timestamp <- format(timing$end_time, "%Y-%m-%d %H:%M:%S")

    # Determine log directory
    if (is.null(timing_log_dir)) {
      if (!is.null(mutantigen_dir)) {
        timing_log_dir <- file.path(mutantigen_dir, "timing_logs")
      } else {
        timing_log_dir <- file.path(here::here(), "mutantigen_parallel", "timing_logs")
      }
    }

    # Create timing log directory if it doesn't exist
    if (!dir.exists(timing_log_dir)) {
      dir.create(timing_log_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Save timing log as RDS (for easy R analysis) and CSV (for human reading)
    timing_file_rds <- file.path(timing_log_dir, paste0(method, "_NUM_", NUM, "_timing.rds"))
    timing_file_csv <- file.path(timing_log_dir, paste0(method, "_NUM_", NUM, "_timing.csv"))

    tryCatch({
      saveRDS(timing, timing_file_rds)

      # Create CSV summary
      csv_data <- data.frame(
        NUM = timing$NUM,
        method = timing$method,
        start_time = timing$start_timestamp,
        end_time = timing$end_timestamp,
        total_time_sec = round(timing$total_time_sec, 2),
        total_time_min = round(timing$total_time_min, 2),
        yaml_time_sec = round(timing$yaml_creation_time_sec %||% NA, 3),
        simulation_time_sec = round(timing$simulation_time_sec %||% NA, 2),
        extraction_time_sec = round(timing$extraction_time_sec %||% NA, 3),
        exit_code = timing$exit_code %||% NA,
        success = timing$success,
        error_message = timing$error_message,
        stringsAsFactors = FALSE
      )

      write.csv(csv_data, timing_file_csv, row.names = FALSE)

      if (verbose) {
        cat("\nTiming log saved to:\n")
        cat("  RDS:", timing_file_rds, "\n")
        cat("  CSV:", timing_file_csv, "\n")
      }
    }, error = function(e) {
      if (verbose) cat("WARNING: Could not save timing log:", conditionMessage(e), "\n")
    })

    # Always restore working directory
    setwd(original_wd)
  })

  return(result)
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a
