# depends: 
analyze_mutantigen_outputs <- function(
  NUM,
  directory = ".",
  peak_col = "totalCases",
  min_peak_height = NULL,
  min_peak_prominence = NULL,
  min_peak_distance_rows = 1,
  waterline_method = c("median"),
  assign_samples_to_peak = c("nearest"),
  peak_window = NULL,
  nonrare_type_threshold = 100
) {
  ##########################################################################
  # analyze_mutantigen_NUM
  #
  # PURPOSE
  #   Read MutAntiGen output files for simulation NUM and compute:
  #
  #   1) Epidemic peak summaries from out_NUM.timeseries using totalCases:
  #      - number of peaks
  #      - times of peaks
  #      - time between peaks
  #      - variation in time between peaks
  #      - peak heights
  #      - variation in peak heights
  #      - peak widths (full width at half maximum, FWHM)
  #      - "water line" where half of totalCases values are below it (median)
  #      - number of peaks above that water line
  #      - total infection burden from totalCases
  #
  #   2) Antigenic-type summaries from sampled viruses:
  #      - extract terminal antigenic type for each sampled tip from
  #        out_NUM.simmapAntigenic
  #      - assign sampled tips to epidemic peaks by time
  #      - determine the dominant antigenic type at each peak
  #      - determine whether dominant type changes from peak to peak
  #      - count number of distinct dominant types across peaks
  #      - count total antigenic types observed among sampled tips
  #      - count non-rare antigenic types among sampled tips
  #
  #   3) Load / fitness-proxy summaries:
  #      - extract terminal mutational load for each sampled tip from
  #        out_NUM.simmapLoad
  #      - summarize load overall and by peak
  #      - return tip-level table with sample time, antigenic type, load
  #
  #   4) Population mutation-load summaries from out_NUM.mutationSeries:
  #      - summarize mean and variance of mutation load through time
  #
  #
  # FILE USAGE
  #
  #   out_NUM.timeseries
  #     Used for the epidemic curve on peak_col = "totalCases" by default.
  #     date = simulation time
  #     totalCases = incidence-like measure saved through time
  #
  #   out_NUM.tips
  #     Used for sampled viruses:
  #       name = sampled tip id
  #       year = sample time
  #       ag1, ag2 = antigenic coordinates
  #
  #   out_NUM.simmapAntigenic
  #     SimMap annotation of antigenic type along branches.
  #     For each tip name found in out_NUM.tips, we extract the FINAL state
  #     on the corresponding annotated branch and treat it as that tip's
  #     antigenic type at sampling time.
  #
  #   out_NUM.simmapLoad
  #     SimMap annotation of mutational load along branches.
  #     For each tip name found in out_NUM.tips, we extract the FINAL state
  #     on the corresponding annotated branch and treat it as that tip's
  #     load at sampling time.
  #
  #   out_NUM.mutationSeries
  #     Each row is interpreted as the distribution of the population across
  #     mutational load classes at a time point. The columns are load classes
  #     0,1,2,... and the values are counts in each class.
  #     This function computes row-wise mean and variance of mutation load.
  #
  #
  # IMPORTANT LIMITATIONS
  #
  #   - "Dominant antigenic type per peak" is SAMPLE-BASED, not full-population.
  #     It is based on sampled tips assigned to each peak by time.
  #
  #   - "Number of antigenic types with >= 100 cases" cannot be computed exactly
  #     from these files unless there is a file with type-specific case counts.
  #     Here, we compute the number of antigenic types with >= threshold SAMPLED
  #     TIPS instead.
  #
  #   - "Absolute fitness for each virus" is not directly available from these
  #     files. We return mutational load from simmapLoad as a fitness-related
  #     proxy. Interpretation depends on the simulation design.
  #
  #   - "Total number of people infected by antigenic type" cannot be computed
  #     exactly from these files alone. We return sampled-tip counts by type.
  #
  ##########################################################################

  waterline_method <- match.arg(waterline_method)
  assign_samples_to_peak <- match.arg(assign_samples_to_peak)

  # ---------------------------
  # Helper: safe path builder
  # ---------------------------
  make_path <- function(suffix) {
    file.path(directory, paste0("out_", NUM, suffix))
  }

  # ---------------------------
  # Helper: local maxima
  # ---------------------------
  local_maxima <- function(y) {
    n <- length(y)
    if (n < 3) return(integer(0))
    which(y[2:(n - 1)] > y[1:(n - 2)] & y[2:(n - 1)] >= y[3:n]) + 1L
  }

  # ---------------------------
  # Helper: approximate prominence
  # ---------------------------
  approx_prominence <- function(y, i) {
    n <- length(y)

    l <- i - 1L
    while (l > 1L && y[l] < y[i]) l <- l - 1L
    left_min <- min(y[l:i], na.rm = TRUE)

    r <- i + 1L
    while (r < n && y[r] <= y[i]) r <- r + 1L
    right_min <- min(y[i:min(r, n)], na.rm = TRUE)

    y[i] - max(left_min, right_min)
  }

  # ---------------------------
  # Helper: filter peaks
  # ---------------------------
  filter_peaks <- function(y, candidates, min_height, min_prominence, min_distance_rows) {
    if (length(candidates) == 0) return(integer(0))

    ord <- candidates[order(y[candidates], decreasing = TRUE)]
    selected <- integer(0)

    for (i in ord) {
      if (!is.null(min_height) && y[i] < min_height) next
      if (!is.null(min_prominence)) {
        prom <- approx_prominence(y, i)
        if (prom < min_prominence) next
      }
      if (length(selected) > 0 && any(abs(selected - i) < min_distance_rows)) next
      selected <- c(selected, i)
    }

    sort(selected)
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
  # Helper: parse a SimMap state string
  #
  # Example:
  #   "{1,0.0767:1,0.0082:1,0.0465:2}"
  #
  # Interpretation used here:
  #   The final state is the value after the last colon,
  #   or the initial value if there are no transitions.
  # ---------------------------
  parse_final_state <- function(state_string) {
    s <- gsub("^\\{|\\}$", "", state_string)
    s <- trimws(s)

    # no colon => single state
    if (!grepl(":", s, fixed = TRUE)) {
      suppressWarnings(val <- as.numeric(s))
      return(val)
    }

    pieces <- strsplit(s, ",", fixed = TRUE)[[1]]
    last_piece <- tail(pieces, 1)
    final_state_str <- sub("^.*:", "", last_piece)
    suppressWarnings(as.numeric(final_state_str))
  }

  # ---------------------------
  # Helper: parse all node/state pairs from SimMap line
  #
  # Extract patterns like:
  #   name:{...}
  #
  # Returns data.frame(name, final_state)
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
    finals <- vapply(states_vec, parse_final_state, numeric(1))

    out <- data.frame(
      name = names_vec,
      final_state = finals,
      stringsAsFactors = FALSE
    )

    # if duplicated names appear, keep the last occurrence
    out <- out[!duplicated(out$name, fromLast = TRUE), ]
    rownames(out) <- NULL
    out
  }

  # ---------------------------
  # Helper: summarize numeric vector
  # ---------------------------
  num_summary <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(list(n = 0, mean = NA_real_, sd = NA_real_, var = NA_real_, cv = NA_real_,
                  min = NA_real_, median = NA_real_, max = NA_real_))
    }
    m <- mean(x)
    s <- if (length(x) >= 2) sd(x) else NA_real_
    list(
      n = length(x),
      mean = m,
      sd = s,
      var = if (length(x) >= 2) var(x) else NA_real_,
      cv = if (!is.na(s) && !isTRUE(all.equal(m, 0))) s / m else NA_real_,
      min = min(x),
      median = median(x),
      max = max(x)
    )
  }

  # ---------------------------
  # 1) Read timeseries
  # ---------------------------
  timeseries_path <- make_path(".timeseries")
  ts <- read.table(timeseries_path, header = TRUE, stringsAsFactors = FALSE)

  if (!("date" %in% names(ts))) stop("timeseries file must contain a 'date' column")
  if (!(peak_col %in% names(ts))) stop(sprintf("timeseries file must contain '%s'", peak_col))

  x <- ts$date
  y <- ts[[peak_col]]

  # ---------------------------
  # 2) Detect peaks in totalCases
  # ---------------------------
  candidate_peaks <- local_maxima(y)
  peaks_idx <- filter_peaks(
    y = y,
    candidates = candidate_peaks,
    min_height = min_peak_height,
    min_prominence = min_peak_prominence,
    min_distance_rows = min_peak_distance_rows
  )

  peak_table <- data.frame(
    peak_index = integer(0),
    peak_time = numeric(0),
    peak_height = numeric(0),
    left_half_height_time = numeric(0),
    right_half_height_time = numeric(0),
    peak_width_fwhm = numeric(0),
    stringsAsFactors = FALSE
  )

  if (length(peaks_idx) > 0) {
    widths <- t(sapply(peaks_idx, function(i) peak_fwhm(x, y, i)))
    peak_table <- data.frame(
      peak_index = peaks_idx,
      peak_time = x[peaks_idx],
      peak_height = y[peaks_idx],
      left_half_height_time = widths[, "left_cross"],
      right_half_height_time = widths[, "right_cross"],
      peak_width_fwhm = widths[, "width"],
      stringsAsFactors = FALSE
    )
  }

  interpeak_times <- if (nrow(peak_table) >= 2) diff(peak_table$peak_time) else numeric(0)

  waterline <- switch(
    waterline_method,
    median = median(y, na.rm = TRUE)
  )

  peaks_above_waterline <- sum(peak_table$peak_height > waterline, na.rm = TRUE)

  epidemic_summary <- list(
    description = paste(
      "Computed from out_", NUM, ".timeseries using column ", peak_col, ". ",
      "Peaks are local maxima filtered by optional minimum height, optional minimum prominence, ",
      "and minimum row distance. Peak width is full width at half maximum (FWHM). ",
      "Water line is the median of ", peak_col, ", i.e. the horizontal level below which half ",
      "of the observed ", peak_col, " values fall.",
      sep = ""
    ),
    number_of_peaks = nrow(peak_table),
    peak_table = peak_table,
    interpeak_times = interpeak_times,
    interpeak_time_summary = num_summary(interpeak_times),
    peak_height_summary = num_summary(peak_table$peak_height),
    peak_width_summary = num_summary(peak_table$peak_width_fwhm),
    waterline = waterline,
    peaks_above_waterline = peaks_above_waterline,
    total_cases_sum = sum(ts[[peak_col]], na.rm = TRUE),
    total_cases_trapz = {
      if (length(x) >= 2) {
        sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
  )

  # ---------------------------
  # 3) Read tips
  # ---------------------------
  tips_path <- make_path(".tips")
  tips <- read.csv(tips_path, stringsAsFactors = FALSE)

  required_tip_cols <- c("name", "year")
  missing_tip_cols <- setdiff(required_tip_cols, names(tips))
  if (length(missing_tip_cols) > 0) {
    stop("tips file is missing required columns: ", paste(missing_tip_cols, collapse = ", "))
  }

  # ---------------------------
  # 4) Parse simmapAntigenic and simmapLoad
  # ---------------------------
  simmap_antigenic_path <- make_path(".simmapAntigenic")
  simmap_load_path <- make_path(".simmapLoad")

  simmap_antigenic <- parse_simmap_file(simmap_antigenic_path)
  names(simmap_antigenic)[names(simmap_antigenic) == "final_state"] <- "antigenic_type"

  simmap_load <- parse_simmap_file(simmap_load_path)
  names(simmap_load)[names(simmap_load) == "final_state"] <- "mutational_load"

  # Join terminal antigenic type and terminal load onto sampled tips
  tip_states <- merge(tips, simmap_antigenic, by = "name", all.x = TRUE)
  tip_states <- merge(tip_states, simmap_load, by = "name", all.x = TRUE)

  # ---------------------------
  # 5) Assign sampled tips to peaks by time
  #
  # Default:
  #   nearest peak in time, optionally restricted by peak_window
  #
  # If peak_window is NULL:
  #   use half the minimum interpeak distance as a default,
  #   or Inf if there is only one peak.
  # ---------------------------
  if (is.null(peak_window)) {
    if (length(interpeak_times) >= 1) {
      peak_window <- min(interpeak_times) / 2
    } else {
      peak_window <- Inf
    }
  }

  assign_tip_to_peak <- function(sample_time, peak_times, peak_window) {
    if (length(peak_times) == 0) return(NA_integer_)
    d <- abs(peak_times - sample_time)
    j <- which.min(d)
    if (d[j] <= peak_window) return(j)
    NA_integer_
  }

  tip_states$peak_id <- if (nrow(peak_table) > 0) {
    vapply(tip_states$year,
           assign_tip_to_peak,
           integer(1),
           peak_times = peak_table$peak_time,
           peak_window = peak_window)
  } else {
    rep(NA_integer_, nrow(tip_states))
  }

  # ---------------------------
  # 6) Antigenic-type summaries
  # ---------------------------

  # Total observed antigenic types among sampled tips
  observed_antigenic_types <- sort(unique(tip_states$antigenic_type[!is.na(tip_states$antigenic_type)]))
  total_observed_antigenic_types <- length(observed_antigenic_types)

  # Non-rare antigenic types among sampled tips
  type_counts <- sort(table(tip_states$antigenic_type), decreasing = TRUE)
  type_counts <- type_counts[names(type_counts) != "NA"]
  nonrare_type_count <- sum(type_counts >= nonrare_type_threshold)

  # Dominant antigenic type for each peak, based on sampled tips assigned to that peak
  peak_antigenic_summary <- data.frame(
    peak_id = integer(0),
    peak_time = numeric(0),
    n_sampled_tips = integer(0),
    dominant_antigenic_type = numeric(0),
    dominant_type_tip_count = integer(0),
    dominant_type_fraction = numeric(0),
    n_types_observed_in_peak_samples = integer(0),
    stringsAsFactors = FALSE
  )

  if (nrow(peak_table) > 0) {
    peak_antigenic_summary <- do.call(
      rbind,
      lapply(seq_len(nrow(peak_table)), function(k) {
        sub <- tip_states[tip_states$peak_id == k & !is.na(tip_states$antigenic_type), , drop = FALSE]
        if (nrow(sub) == 0) {
          data.frame(
            peak_id = k,
            peak_time = peak_table$peak_time[k],
            n_sampled_tips = 0L,
            dominant_antigenic_type = NA_real_,
            dominant_type_tip_count = 0L,
            dominant_type_fraction = NA_real_,
            n_types_observed_in_peak_samples = 0L,
            stringsAsFactors = FALSE
          )
        } else {
          tab <- sort(table(sub$antigenic_type), decreasing = TRUE)
          dom_type <- as.numeric(names(tab)[1])
          dom_n <- as.integer(tab[1])
          data.frame(
            peak_id = k,
            peak_time = peak_table$peak_time[k],
            n_sampled_tips = nrow(sub),
            dominant_antigenic_type = dom_type,
            dominant_type_tip_count = dom_n,
            dominant_type_fraction = dom_n / nrow(sub),
            n_types_observed_in_peak_samples = length(tab),
            stringsAsFactors = FALSE
          )
        }
      })
    )
  }

  dominant_types_across_peaks <- peak_antigenic_summary$dominant_antigenic_type
  dominant_types_across_peaks <- dominant_types_across_peaks[!is.na(dominant_types_across_peaks)]

  antigenic_summary <- list(
    description = paste(
      "Antigenic summaries are sample-based and use three files together: ",
      "out_", NUM, ".tips provides sampled virus names and sample times; ",
      "out_", NUM, ".simmapAntigenic provides branch annotations of antigenic type, from which ",
      "the final state on each sampled tip branch is used as that tip's antigenic type; ",
      "tips are assigned to the nearest epidemic peak in time using peak times from out_", NUM,
      ".timeseries. Dominant antigenic type per peak means the most frequent antigenic type among ",
      "sampled tips assigned to that peak. Counts of antigenic types here refer to sampled tips, ",
      "not full-population case counts.",
      sep = ""
    ),
    total_observed_antigenic_types_in_sampled_tips = total_observed_antigenic_types,
    observed_antigenic_types = observed_antigenic_types,
    antigenic_type_tip_counts = type_counts,
    nonrare_antigenic_type_threshold_in_sampled_tips = nonrare_type_threshold,
    nonrare_antigenic_type_count_in_sampled_tips = nonrare_type_count,
    peak_antigenic_summary = peak_antigenic_summary,
    number_of_distinct_dominant_types_across_peaks = length(unique(dominant_types_across_peaks)),
    dominant_type_changes_between_adjacent_peaks = {
      if (length(dominant_types_across_peaks) >= 2) {
        sum(diff(dominant_types_across_peaks) != 0, na.rm = TRUE)
      } else {
        0L
      }
    }
  )

  # ---------------------------
  # 7) Load / fitness-proxy summaries
  # ---------------------------
  peak_load_summary <- data.frame(
    peak_id = integer(0),
    peak_time = numeric(0),
    n_sampled_tips_with_load = integer(0),
    mean_load = numeric(0),
    var_load = numeric(0),
    sd_load = numeric(0),
    stringsAsFactors = FALSE
  )

  if (nrow(peak_table) > 0) {
    peak_load_summary <- do.call(
      rbind,
      lapply(seq_len(nrow(peak_table)), function(k) {
        sub <- tip_states[tip_states$peak_id == k & !is.na(tip_states$mutational_load), , drop = FALSE]
        if (nrow(sub) == 0) {
          data.frame(
            peak_id = k,
            peak_time = peak_table$peak_time[k],
            n_sampled_tips_with_load = 0L,
            mean_load = NA_real_,
            var_load = NA_real_,
            sd_load = NA_real_,
            stringsAsFactors = FALSE
          )
        } else {
          data.frame(
            peak_id = k,
            peak_time = peak_table$peak_time[k],
            n_sampled_tips_with_load = nrow(sub),
            mean_load = mean(sub$mutational_load, na.rm = TRUE),
            var_load = if (nrow(sub) >= 2) var(sub$mutational_load, na.rm = TRUE) else NA_real_,
            sd_load = if (nrow(sub) >= 2) sd(sub$mutational_load, na.rm = TRUE) else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      })
    )
  }

  load_summary <- list(
    description = paste(
      "Load summaries use out_", NUM, ".tips and out_", NUM, ".simmapLoad. ",
      "The final state on each sampled tip branch in simmapLoad is used as that tip's mutational load. ",
      "This is returned as a per-tip value and summarized overall and by epidemic peak. ",
      "It should be interpreted as a mutational-load summary or a fitness-related proxy, not necessarily ",
      "the exact absolute fitness variable unless confirmed from the simulation code.",
      sep = ""
    ),
    overall_tip_load_summary = num_summary(tip_states$mutational_load),
    peak_load_summary = peak_load_summary
  )

  # ---------------------------
  # 8) mutationSeries summaries
  # ---------------------------
  mutation_series_path <- make_path(".mutationSeries")
  mutation_series_summary <- NULL

  if (file.exists(mutation_series_path)) {
    ms <- read.table(mutation_series_path, sep = ",", header = FALSE, strip.white = TRUE)
    load_classes <- 0:(ncol(ms) - 1)

    row_mean_load <- apply(ms, 1, function(counts) {
      counts <- as.numeric(counts)
      total <- sum(counts)
      if (total == 0) return(NA_real_)
      sum(load_classes * counts) / total
    })

    row_var_load <- apply(ms, 1, function(counts) {
      counts <- as.numeric(counts)
      total <- sum(counts)
      if (total == 0) return(NA_real_)
      mu <- sum(load_classes * counts) / total
      sum(((load_classes - mu)^2) * counts) / total
    })

    mutation_series_summary <- list(
      description = paste(
        "Mutation-series summaries use out_", NUM, ".mutationSeries. ",
        "Each row is treated as the population distribution across mutational load classes 0,1,2,... . ",
        "For each row, the function computes the population mean and variance of mutational load. ",
        "These are population-level summaries over time, not tip-specific values.",
        sep = ""
      ),
      n_timepoints = nrow(ms),
      n_load_classes = ncol(ms),
      row_mean_load = row_mean_load,
      row_var_load = row_var_load,
      overall_mean_of_row_means = mean(row_mean_load, na.rm = TRUE),
      overall_mean_of_row_vars = mean(row_var_load, na.rm = TRUE)
    )
  }

  # ---------------------------
  # 9) Sample-based type burden summaries
  # ---------------------------
  # This is NOT true infections by type. It is sampled tips by type.
  sampled_type_burden_summary <- NULL
  if (length(type_counts) > 0) {
    tc <- as.numeric(type_counts)
    sampled_type_burden_summary <- list(
      description = paste(
        "This section summarizes the distribution of sampled tips across antigenic types. ",
        "It uses out_", NUM, ".tips joined to terminal antigenic states from out_", NUM,
        ".simmapAntigenic. These are sampled-tip counts by type, not true counts of all infected hosts.",
        sep = ""
      ),
      n_types = length(tc),
      total_sampled_tips_with_antigenic_type = sum(tc),
      type_count_summary = num_summary(tc),
      max_minus_min_type_count = max(tc) - min(tc),
      top_type_fraction = max(tc) / sum(tc)
    )
  }

  # ---------------------------
  # Final return object
  # ---------------------------
  result <- list(
    input_files = list(
      timeseries = timeseries_path,
      tips = tips_path,
      simmapAntigenic = simmap_antigenic_path,
      simmapLoad = simmap_load_path,
      mutationSeries = if (file.exists(mutation_series_path)) mutation_series_path else NA_character_
    ),
    interpretation_notes = c(
      "Peak summaries are computed on totalCases by default.",
      "Peak widths are full width at half maximum (FWHM).",
      "The water line is the median totalCases value through time.",
      "Antigenic dominance by peak is based on sampled tips assigned to the nearest epidemic peak by time.",
      "Non-rare antigenic types are counted using sampled tips, not cases.",
      "Mutational load from simmapLoad is returned as a fitness-related proxy, not guaranteed to be exact absolute fitness.",
      "True total infected by antigenic type cannot be recovered exactly from these files alone."
    ),
    epidemic_summary = epidemic_summary,
    antigenic_summary = antigenic_summary,
    load_summary = load_summary,
    mutation_series_summary = mutation_series_summary,
    sampled_type_burden_summary = sampled_type_burden_summary,
    tip_level_data = tip_states
  )

  class(result) <- "mutantigen_analysis"
  return(result)
}