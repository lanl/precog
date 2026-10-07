# depends: 
extract_mutantigen_10_outputs <- function(
  NUM,
  directory = ".",
  first_month_duration = 0.077,  # ~4 weeks in years (28/365), to exclude at the beginning
  peak_threshold_pct = 0.05,      # peaks must exceed 5% of population
  type_threshold_pct = 0.01       # types must have >=1% of population in samples
) {
  ##########################################################################
  # extract_mutantigen_10_outputs
  #
  # PURPOSE
  #   Extract 10 specific summary statistics from MutAntiGen
  #   simulation outputs for use in output-space expansion.
  #
  # OUTPUTS (in order):
  #   1. Total number of antigenic types with >= 1% of population size in sampled tips
  #   2. Number of unique antigens with dominant proportions across all peaks (within FWHM)
  #   3. Number of peaks > 5% of total population size
  #   4. Average meanR across viralFitnessSeries (excluding first month)
  #   5. Average meanR over FWHM regions of peaks only
  #   6. Number of times dominant antigenic type changes from peak to peak
  #   7. Max height of a peak
  #   8. Average number of unique antigens in FWHM region of a peak
  #   9. Total infection burden
  #  10. Range of peak heights (max - min)
  #
  # NOTES
  #   - All calculations exclude the first month of data
  #   - FWHM = Full Width at Half Maximum
  #   - Peaks are filtered to only include those > 5% of population size
  #   - Population size is read from totalN column in .timeseries
  #   - R values come from meanR in .viralFitnessSeries
  ##########################################################################

  # ---------------------------
  # Helper: safe path builder
  # ---------------------------
  make_path <- function(suffix) {
    file.path(directory, paste0("out_", NUM, suffix))
  }

  # ---------------------------
  # Read timeseries and filter out first month
  # ---------------------------
  timeseries_path <- make_path(".timeseries")
  ts_raw <- read.table(timeseries_path, header = TRUE, stringsAsFactors = FALSE)

  if (!("date" %in% names(ts_raw))) stop("timeseries file must contain a 'date' column")
  if (!("totalCases" %in% names(ts_raw))) stop("timeseries file must contain 'totalCases'")
  if (!("totalN" %in% names(ts_raw))) stop("timeseries file must contain 'totalN'")

  # Extract population size (assume constant or use median)
  population_size <- median(ts_raw$totalN, na.rm = TRUE)

  # Filter out first month
  ts <- ts_raw[ts_raw$date > first_month_duration, ]

  if (nrow(ts) < 3) {
    warning("Not enough data after filtering first month")
    return(rep(NA_real_, 10))
  }

  x <- ts$date
  y <- ts$totalCases

  # ---------------------------
  # Read viralFitnessSeries and filter out first month
  # ---------------------------
  viral_fitness_path <- make_path(".viralFitnessSeries")
  vf_raw <- read.table(viral_fitness_path, header = TRUE, stringsAsFactors = FALSE)

  if (!("date" %in% names(vf_raw))) stop("viralFitnessSeries file must contain a 'date' column")
  if (!("meanR" %in% names(vf_raw))) stop("viralFitnessSeries file must contain 'meanR'")

  # Filter out first month
  vf <- vf_raw[vf_raw$date > first_month_duration, ]

  if (nrow(vf) < 1) {
    warning("No viralFitnessSeries data after filtering first month")
    return(rep(NA_real_, 10))
  }

  x_vf <- vf$date
  y_vf <- vf$meanR

  # ---------------------------
  # Helper: local maxima
  # ---------------------------
  local_maxima <- function(y) {
    n <- length(y)
    if (n < 3) return(integer(0))
    which(y[2:(n - 1)] > y[1:(n - 2)] & y[2:(n - 1)] >= y[3:n]) + 1L
  }

  # ---------------------------
  # Helper: linear interpolation
  # ---------------------------
  interp_x <- function(x0, y0, x1, y1, target) {
    if (isTRUE(all.equal(y0, y1))) return((x0 + x1) / 2)
    x0 + (target - y0) * (x1 - x0) / (y1 - y0)
  }

  # ---------------------------
  # Helper: FWHM for one peak
  # ---------------------------
  peak_fwhm <- function(x, y, i) {
    half <- y[i] / 2

    l <- i
    while (l > 1L && y[l] >= half) l <- l - 1L
    if (l == i) {
      left_cross <- x[i]
    } else {
      left_cross <- interp_x(x[l], y[l], x[l + 1L], y[l + 1L], half)
    }

    r <- i
    while (r < length(y) && y[r] >= half) r <- r + 1L
    if (r == i) {
      right_cross <- x[i]
    } else if (r > length(y)) {
      right_cross <- x[length(y)]
    } else {
      right_cross <- interp_x(x[r - 1L], y[r - 1L], x[r], y[r], half)
    }

    c(left_cross = left_cross,
      right_cross = right_cross,
      width = right_cross - left_cross)
  }

  # ---------------------------
  # Helper: parse SimMap final state (VECTORIZED)
  # ---------------------------
  # Process all state strings at once instead of one-by-one
  # This is 5-10x faster for large files because:
  #   - gsub/sub/grepl operate on vectors in compiled C code
  #   - No R function call overhead per string
  #   - No list creation from strsplit
  parse_final_state_vectorized <- function(state_strings) {
    # Remove { and } from all strings at once
    cleaned <- gsub("^\\{|\\}$", "", state_strings)
    cleaned <- trimws(cleaned)

    # Check which strings have transitions (contain ":")
    has_colon <- grepl(":", cleaned, fixed = TRUE)

    # Initialize result vector
    n <- length(state_strings)
    result <- numeric(n)

    # For strings without colons: entire string is the state
    if (any(!has_colon)) {
      result[!has_colon] <- suppressWarnings(as.numeric(cleaned[!has_colon]))
    }

    # For strings with colons: extract final state after last comma or colon
    # Regex: ".*[,:]([^,:]+)$" captures everything after the last delimiter
    if (any(has_colon)) {
      finals <- sub(".*[,:]([^,:]+)$", "\\1", cleaned[has_colon])
      result[has_colon] <- suppressWarnings(as.numeric(finals))
    }

    return(result)
  }

  # ---------------------------
  # Helper: parse SimMap file
  # ---------------------------
  parse_simmap_file <- function(path) {
    txt <- paste(readLines(path, warn = FALSE), collapse = "")
    matches <- gregexpr("([A-Za-z0-9]+):\\{[^\\}]+\\}", txt, perl = TRUE)[[1]]

    if (length(matches) == 1L && matches[1] == -1L) {
      return(data.frame(name = character(0), final_state = numeric(0), stringsAsFactors = FALSE))
    }

    hits <- regmatches(txt, gregexpr("([A-Za-z0-9]+):\\{[^\\}]+\\}", txt, perl = TRUE))[[1]]

    names_vec <- sub(":\\{.*$", "", hits, perl = TRUE)
    states_vec <- sub("^[A-Za-z0-9]+:(\\{.*\\})$", "\\1", hits, perl = TRUE)

    # Use vectorized parsing instead of vapply loop
    finals <- parse_final_state_vectorized(states_vec)

    out <- data.frame(
      name = names_vec,
      final_state = finals,
      stringsAsFactors = FALSE
    )

    out <- out[!duplicated(out$name, fromLast = TRUE), ]
    rownames(out) <- NULL
    out
  }

  # ---------------------------
  # Detect peaks (>5% of population size)
  # ---------------------------
  min_peak_height <- peak_threshold_pct * population_size
  candidate_peaks <- local_maxima(y)
  peaks_idx <- candidate_peaks[y[candidate_peaks] > min_peak_height]

  if (length(peaks_idx) == 0) {
    warning("No peaks found above threshold")
    return(rep(NA_real_, 10))
  }

  # Calculate FWHM for each peak
  peak_data <- lapply(peaks_idx, function(i) {
    fwhm <- peak_fwhm(x, y, i)
    list(
      idx = i,
      time = x[i],
      height = y[i],
      left = fwhm["left_cross"],
      right = fwhm["right_cross"],
      width = fwhm["width"]
    )
  })

  peak_times <- sapply(peak_data, function(p) p$time)
  peak_heights <- sapply(peak_data, function(p) p$height)

  # ---------------------------
  # Output 3: Number of peaks
  # ---------------------------
  out_3 <- length(peaks_idx)

  # ---------------------------
  # Output 7: Max peak height
  # ---------------------------
  out_7 <- max(peak_heights)

  # ---------------------------
  # Output 10: Range of peak heights
  # ---------------------------
  out_10_range <- max(peak_heights) - min(peak_heights)

  # ---------------------------
  # Output 4: Average meanR across viralFitnessSeries (minus first month)
  # ---------------------------
  out_4 <- mean(y_vf, na.rm = TRUE)

  # ---------------------------
  # Output 5: Average meanR over FWHM regions only
  # ---------------------------
  # For each viralFitnessSeries timepoint, check if it falls in any FWHM region
  fwhm_mask_vf <- rep(FALSE, length(x_vf))
  for (p in peak_data) {
    fwhm_mask_vf <- fwhm_mask_vf | (x_vf >= p$left & x_vf <= p$right)
  }
  out_5 <- if (any(fwhm_mask_vf)) {
    mean(y_vf[fwhm_mask_vf], na.rm = TRUE)
  } else {
    NA_real_
  }

  # ---------------------------
  # Output 9: Total infection burden
  # ---------------------------
  out_9_burden <- if (length(x) >= 2) {
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
  } else {
    sum(y, na.rm = TRUE)
  }

  # ---------------------------
  # Read tips and antigenic data
  # ---------------------------
  tips_path <- make_path(".tips")
  tips <- read.csv(tips_path, stringsAsFactors = FALSE)

  # Filter tips to only those after first month
  tips <- tips[tips$year > first_month_duration, ]

  if (nrow(tips) == 0) {
    warning("No tips found after filtering first month")
    return(c(NA_real_, NA_real_, out_3, out_4, out_5,
             NA_real_, out_7, NA_real_, out_9_burden, out_10_range))
  }

  simmap_antigenic_path <- make_path(".out.simmapAntigenic")
  simmap_antigenic <- parse_simmap_file(simmap_antigenic_path)
  names(simmap_antigenic)[names(simmap_antigenic) == "final_state"] <- "antigenic_type"

  tip_states <- merge(tips, simmap_antigenic, by = "name", all.x = TRUE)
  tip_states <- tip_states[!is.na(tip_states$antigenic_type), ]

  # ---------------------------
  # Output 1: Number of antigenic types >= 1% of population
  # ---------------------------
  min_type_count <- type_threshold_pct * population_size
  type_counts <- table(tip_states$antigenic_type)
  out_1 <- sum(type_counts >= min_type_count)

  # ---------------------------
  # Assign tips to peaks (within FWHM regions)
  # ---------------------------
  tip_states$peak_id <- NA_integer_
  tip_states$in_fwhm <- FALSE

  for (k in seq_along(peak_data)) {
    p <- peak_data[[k]]
    in_region <- tip_states$year >= p$left & tip_states$year <= p$right
    tip_states$peak_id[in_region] <- k
    tip_states$in_fwhm[in_region] <- TRUE
  }

  # ---------------------------
  # Calculate dominant type per peak (within FWHM)
  # ---------------------------
  dominant_types <- integer(0)
  unique_types_per_peak <- integer(0)

  for (k in seq_along(peak_data)) {
    peak_tips <- tip_states[tip_states$peak_id == k & tip_states$in_fwhm, ]

    if (nrow(peak_tips) > 0) {
      type_tab <- table(peak_tips$antigenic_type)
      type_props <- type_tab / sum(type_tab)

      # Dominant type (highest proportion)
      dom_type <- as.numeric(names(which.max(type_props)))
      dominant_types <- c(dominant_types, dom_type)

      # Number of unique types in this peak's FWHM
      unique_types_per_peak <- c(unique_types_per_peak, length(type_tab))
    } else {
      dominant_types <- c(dominant_types, NA_integer_)
      unique_types_per_peak <- c(unique_types_per_peak, 0L)
    }
  }

  # ---------------------------
  # Output 2: Number of unique antigens with dominant proportions across peaks
  # ---------------------------
  out_2 <- length(unique(dominant_types[!is.na(dominant_types)]))

  # ---------------------------
  # Output 6: Number of times dominant type changes
  # ---------------------------
  dom_valid <- dominant_types[!is.na(dominant_types)]
  out_6 <- if (length(dom_valid) >= 2) {
    sum(diff(dom_valid) != 0)
  } else {
    0L
  }

  # ---------------------------
  # Output 8: Average number of unique antigens in FWHM regions
  # ---------------------------
  out_8_unique <- mean(unique_types_per_peak)

  # ---------------------------
  # Return vector of 10 outputs
  # ---------------------------
  return(c(
    out_1,         # 1. Total antigenic types >= 1% pop
    out_2,         # 2. Unique dominant types across peaks
    out_3,         # 3. Number of peaks
    out_4,         # 4. Average meanR (timeseries)
    out_5,         # 5. Average meanR (FWHM only)
    out_6,         # 6. Dominant type changes
    out_7,         # 7. Max peak height
    out_8_unique,  # 8. Average unique antigens per peak FWHM
    out_9_burden,  # 9. Total infection burden
    out_10_range   # 10. Range of peak heights
  ))
}
