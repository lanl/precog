# depends: 
#' @export
assign_regime2 <- function(trial, snippet_length, global_slope_thres) {
  vec <- trial
  regime <- return_best_shapelet_pearson(vec, snippet_length, global_slope_thres)
  return(regime)
}

#' @export
similarity_matrix <- function(vector1, vector2) {
  # Handle near-constant vectors to avoid correlation calculation issues
  if (sd(vector1) < 1e-100 || sd(vector2) < 1e-100) {
    similarity_value <- 0
  } else {
    similarity_value <- cor(vector1, vector2, method = 'pearson')
  }
  similarity_value
}

#' @export
return_best_shapelet_pearson <- function(vector_input, snippet_length, global_slope_thres) {
  # Define names for each shapelet pattern
  shapelet_standard_names <- c("Inc", "Dec", "Surge", "Near Peak", "Crash")
  corrs <- return_all_shapelet_pearson(vector_input, snippet_length, global_slope_thres)
  scenario <- which.max(corrs)
  shapelet_standard_names[scenario]
}

#' @export
return_all_shapelet_pearson <- function(vector_input, snippet_length, global_slope_thres) {
  vector_length <- c(0, snippet_length)   # Look 4 weeks ahead in future while defining shapelet
  shapelet_length <- vector_length[1] + vector_length[2]

  # Define names for each shapelet pattern
  shapelet_standard_names <- c("Inc", "Dec", "Surge", "Near Peak", "Crash")

  # Define shapelet patterns directly as a list
  shapelet_standard_array <- list(
    c(1:snippet_length),              # 1: Increasing
    c(snippet_length:1),              # 2: Decreasing
    c((1:snippet_length)^2),          # 3: Surge
    c(-(1:snippet_length)^(-2)),      # 4: Near peak
    c(-(1:snippet_length)^2)          # 5: Crash/Past peak
  )

  corrs <- numeric(length(shapelet_standard_array))

  # Calculate correlation with each shapelet
  for (i in seq_along(shapelet_standard_array)) {
    corrs[i] <- similarity_matrix(shapelet_standard_array[[i]], vector_input)
  }
  corrs
}

#' @export
do_smooth <- function(truth, x, span = 0.25) { 
  # Use larger span for shorter time series to avoid overfitting
  if(length(truth) < 40) {
    span <- 0.55
  }
  lo_tmp <- loess(truth ~ x, span = span)
  smooth_cases <- predict(lo_tmp, x)
  # Ensure non-negative values (cases can't be negative)
  smooth_cases[smooth_cases < 0] <- 0 
  return(smooth_cases)
}

#' @export
classify_entire_timeseries <- function(sir_vector, snippet_length, match_length = 5) {
  # snippet_length: total length to save (e.g., 9 = 5 observed + 4 future)
  # match_length: length used for shapelet matching (e.g., 5)

  # sir_vector = do_smooth(sir_vector, 1:length(sir_vector))

  all_classifies = NULL
  global_dx <- abs(diff(sir_vector))
  global_slope_thres <- quantile(global_dx, 0.2)

  for(i in 1:(length(sir_vector) - snippet_length)){

    tmp_vec = sir_vector[i:(i+snippet_length-1)]
    # Use only the first match_length values for classification
    tmp_vec_for_matching = tmp_vec[1:match_length]
    tmp_assign = assign_regime2(tmp_vec_for_matching, match_length, global_slope_thres)
    if(length(tmp_assign) == 0) tmp_assign = "None"
    all_classifies = rbind(all_classifies,
            data.frame(value = tmp_vec[1], regime = tmp_assign)
    )
  }

  all_classifies = rbind(all_classifies,
            data.frame(value = tmp_vec[(2):length(tmp_vec)], regime = tmp_assign)
    )

  return(all_classifies)
}


# N          = 10000
# s0         = 0.95
# i0         = 0.0003
# r0         = 1 - i0 - s0

# dave_beta  = 10
# dave_gamma = 0.5
# df2         = sir_sirstuff(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
