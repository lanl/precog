library(tidyverse)
library(data.table)

dir <- "./precog/reporting-delay/"
results_dir <- paste0(dir, "results/")

# Input files -------------------------------------------------------------
variant_file <- "emerging_variant_ts.csv"
metric_file <- "all_metric_results.csv"

variant_raw <- fread(paste0(results_dir, variant_file)) %>%
  transmute(
    Location,
    Date = as.Date(Date),
    delay_days,
    pango,
    p_val
  )

metric_results <- fread(paste0(results_dir, metric_file)) %>%
  mutate(Date = as.Date(Date))

# The validation composition should be identical across delay_days. Confirm
# this before reducing the four repeated copies to one value per lineage/day.
composition_check <- variant_raw %>%
  group_by(Location, Date, pango) %>%
  summarise(
    n_values = n_distinct(p_val),
    .groups = "drop"
  ) %>%
  filter(n_values > 1)

if (nrow(composition_check) > 0) {
  stop(
    "At least one Location-Date-pango combination has different p_val values ",
    "across delay_days; inspect `composition_check` before continuing."
  )
}

variant_ts <- variant_raw %>%
  distinct(Location, Date, pango, p_val)

# Resolve the dominant lineage sequentially within each location. On tied days:
# 1. retain the preceding selected dominant if it is among the tied candidates;
# 2. otherwise select the candidate with the largest proportion on the preceding
#    available date;
# 3. use alphabetical order only if a tie remains.
resolve_location <- function(location_data) {
  location_data <- location_data %>% arrange(Date, pango)
  dates <- sort(unique(location_data$Date))

  dominant_pango <- character(length(dates))
  dominant_p_val <- numeric(length(dates))
  delta_p_val <- rep(NA_real_, length(dates))

  for (i in seq_along(dates)) {
    today <- location_data %>% filter(Date == dates[i])
    candidates <- today %>%
      filter(p_val == max(p_val, na.rm = TRUE)) %>%
      arrange(pango)

    if (i > 1 && dominant_pango[i - 1] %in% candidates$pango) {
      selected_pango <- dominant_pango[i - 1]
    } else if (i > 1) {
      yesterday <- location_data %>%
        filter(Date == dates[i - 1], pango %in% candidates$pango) %>%
        select(pango, previous_p_val = p_val)

      selected_pango <- candidates %>%
        select(pango) %>%
        left_join(yesterday, by = "pango") %>%
        mutate(previous_p_val = coalesce(previous_p_val, 0)) %>%
        arrange(desc(previous_p_val), pango) %>%
        slice(1) %>%
        pull(pango)
    } else {
      selected_pango <- candidates$pango[1]
    }

    dominant_pango[i] <- selected_pango
    dominant_p_val[i] <- today %>%
      filter(pango == selected_pango) %>%
      pull(p_val) %>%
      first()

    if (i > 1) {
      previous_values <- location_data %>%
        filter(Date == dates[i - 1], pango == selected_pango) %>%
        pull(p_val)

      previous_p_val <- if (length(previous_values) == 0) 0 else previous_values[1]
      elapsed_days <- as.numeric(dates[i] - dates[i - 1])
      delta_p_val[i] <- (dominant_p_val[i] - previous_p_val) / elapsed_days
    }
  }

  tibble(
    Date = dates,
    pango_dom = dominant_pango,
    p_val = dominant_p_val,
    delta_p_val = delta_p_val
  )
}

dominant_variant_ts <- variant_ts %>%
  group_by(Location) %>%
  group_modify(~ resolve_location(.x)) %>%
  ungroup()

# Flag calendar-date windows around each change in the selected dominant
# variant. The change date is the first "after" day.
add_turnover_flags <- function(location_data) {
  location_data <- location_data %>% arrange(Date)

  change_dates <- location_data %>%
    mutate(dominant_changed = pango_dom != lag(pango_dom)) %>%
    filter(coalesce(dominant_changed, FALSE)) %>%
    pull(Date)

  if (length(change_dates) == 0) {
    return(
      location_data %>%
        mutate(
          turnover_2wk_after = FALSE,
          turnover_1wk_before_after = FALSE,
          turnover_2wk_before_after = FALSE
        )
    )
  }

  location_data %>%
    mutate(
      turnover_2wk_after = map_lgl(
        Date,
        ~ any(.x >= change_dates & .x <= change_dates + 13)
      ),
      turnover_1wk_before_after = map_lgl(
        Date,
        ~ any(.x >= change_dates - 7 & .x <= change_dates + 6)
      ),
      turnover_2wk_before_after = map_lgl(
        Date,
        ~ any(.x >= change_dates - 14 & .x <= change_dates + 13)
      )
    )
}

dominant_variant_ts <- dominant_variant_ts %>%
  group_by(Location) %>%
  group_modify(~ add_turnover_flags(.x)) %>%
  ungroup()

write_csv(
  dominant_variant_ts,
  paste0(results_dir, "dominant_variant_rate_of_change.csv"),
  na = ""
)

# Merge each validation-composition result onto each corresponding delay-day
# metric record.
dominant_variant_w <- metric_results %>%
  left_join(dominant_variant_ts, by = c("Location", "Date")) %>%
  filter(
    Date < as.Date("2023-01-01"),
    Date >= as.Date("2020-11-01")
  ) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant")

write_csv(
  dominant_variant_w,
  paste0(results_dir, "dominant_variant_w_merged.csv"),
  na = ""
)

