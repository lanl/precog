# Run SIR Parts 2, 3, and 4 for multiple sample sizes
# This script runs Basic LHS (Part 2), OSFD Basic (Part 3), and Wang Original OSFD (Part 4)
# for a sequence of sample sizes
# Author: AC Murph
# Date: Aug 2026, updated Sep 2026
# NOTE: Parts 2 and 3 are currently COMMENTED OUT to run only Part 4 (Wang Original)

library(here)
setwd(here::here())

# Define sample sizes to test
NSAMPLES <- seq(1000, 20000, 1000)
cat(sprintf("Running ONLY Part 4 (Wang Original OSFD) for %d different sample sizes\n", length(NSAMPLES)))
cat(sprintf("Sample sizes: %s\n\n", paste(NSAMPLES, collapse = ", ")))

# Track timing and results
results_summary <- data.frame(
  nsuccesses = integer(),
  part = character(),
  status = character(),
  duration_minutes = numeric(),
  stringsAsFactors = FALSE
)

# Main loop over sample sizes
for (n in NSAMPLES) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Processing nsuccesses = %d\n", n))
  cat(sprintf("========================================\n\n"))

  # # --- Part 2: Basic LHS --- [COMMENTED OUT]
  # cat(sprintf("[%s] Starting Part 2 (Basic LHS) with n=%d...\n", Sys.time(), n))
  # start_time <- Sys.time()
  #
  # ret_code_2 <- system2(
  #   command = "Rscript",
  #   args = c(
  #     here::here("R", "run_sir_part2_basic_lhs.R"),
  #     as.character(n)
  #   ),
  #   stdout = here::here("logfiles", sprintf("part2_n%d.Rout", n)),
  #   stderr = here::here("logfiles", sprintf("part2_n%d.Rout", n))
  # )
  #
  # end_time <- Sys.time()
  # duration_2 <- as.numeric(difftime(end_time, start_time, units = "mins"))
  #
  # status_2 <- if (ret_code_2 == 0) "SUCCESS" else "FAILED"
  # cat(sprintf("[%s] Part 2 completed with status: %s (%.2f minutes)\n",
  #             Sys.time(), status_2, duration_2))
  #
  # results_summary <- rbind(results_summary, data.frame(
  #   nsuccesses = n,
  #   part = "Part2_BasicLHS",
  #   status = status_2,
  #   duration_minutes = duration_2,
  #   stringsAsFactors = FALSE
  # ))

  # # --- Part 3: OSFD Basic --- [COMMENTED OUT]
  # cat(sprintf("[%s] Starting Part 3 (OSFD Basic) with n=%d...\n", Sys.time(), n))
  # start_time <- Sys.time()
  #
  # ret_code_3 <- system2(
  #   command = "Rscript",
  #   args = c(
  #     here::here("R", "run_sir_part3_osfd_basic.R"),
  #     as.character(n)
  #   ),
  #   stdout = here::here("logfiles", sprintf("part3_n%d.Rout", n)),
  #   stderr = here::here("logfiles", sprintf("part3_n%d.Rout", n))
  # )
  #
  # end_time <- Sys.time()
  # duration_3 <- as.numeric(difftime(end_time, start_time, units = "mins"))
  #
  # status_3 <- if (ret_code_3 == 0) "SUCCESS" else "FAILED"
  # cat(sprintf("[%s] Part 3 completed with status: %s (%.2f minutes)\n",
  #             Sys.time(), status_3, duration_3))
  #
  # results_summary <- rbind(results_summary, data.frame(
  #   nsuccesses = n,
  #   part = "Part3_OSFD_Basic",
  #   status = status_3,
  #   duration_minutes = duration_3,
  #   stringsAsFactors = FALSE
  # ))

  # --- Part 4: Wang Original OSFD ---
  cat(sprintf("[%s] Starting Part 4 (Wang Original OSFD) with n=%d...\n", Sys.time(), n))
  start_time <- Sys.time()

  ret_code_4 <- system2(
    command = "Rscript",
    args = c(
      here::here("R", "run_sir_WangOrig_osfd_basic.R"),
      as.character(n)
    ),
    stdout = here::here("logfiles", sprintf("part4_n%d.Rout", n)),
    stderr = here::here("logfiles", sprintf("part4_n%d.Rout", n))
  )

  end_time <- Sys.time()
  duration_4 <- as.numeric(difftime(end_time, start_time, units = "mins"))

  status_4 <- if (ret_code_4 == 0) "SUCCESS" else "FAILED"
  cat(sprintf("[%s] Part 4 completed with status: %s (%.2f minutes)\n",
              Sys.time(), status_4, duration_4))

  results_summary <- rbind(results_summary, data.frame(
    nsuccesses = n,
    part = "Part4_WangOrig",
    status = status_4,
    duration_minutes = duration_4,
    stringsAsFactors = FALSE
  ))

  # Report progress
  cat(sprintf("\nCompleted n=%d: Part4=%s (%.1f min)\n",
              n, status_4, duration_4))
}

# Save results summary
cat("\n========================================\n")
cat("All runs complete!\n")
cat("========================================\n\n")

save(results_summary, file = here::here("data", "sir_parts2_3_4_timing_summary.RData"))
write.csv(results_summary, file = here::here("data", "sir_parts2_3_4_timing_summary.csv"),
          row.names = FALSE)

# Print summary table
cat("\nSummary of all runs:\n")
print(results_summary)

# Calculate totals
total_time <- sum(results_summary$duration_minutes)
n_success <- sum(results_summary$status == "SUCCESS")
n_total <- nrow(results_summary)

cat(sprintf("\nTotal time: %.2f minutes (%.2f hours)\n", total_time, total_time / 60))
cat(sprintf("Success rate: %d/%d (%.1f%%)\n", n_success, n_total, 100 * n_success / n_total))

if (any(results_summary$status == "FAILED")) {
  cat("\nFailed runs:\n")
  print(results_summary[results_summary$status == "FAILED", ])
}

cat("\nResults summary saved to:\n")
cat(sprintf("  - data/sir_parts2_3_4_timing_summary.RData\n"))
cat(sprintf("  - data/sir_parts2_3_4_timing_summary.csv\n"))
cat("\nIndividual log files saved to logfiles/part4_n{N}.Rout\n")
cat("\nNOTE: Parts 2 and 3 were SKIPPED (commented out) - only Part 4 ran.\n")
