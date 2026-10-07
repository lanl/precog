# Run SIR Parts 1, 2, 3, and 4 for a single sample size
# Designed to be called from a SLURM array job
# Usage: Rscript run_sir_parts2_3_single_n.R <nsuccesses>
# Part 1: SIR Maps (inverse PIV/PIT mapping)
# Part 2: Basic LHS (input-space filling)
# Parts 3 and 4: OSFD (output-space filling, nsuccesses used as max_budget)
# Author: AC Murph
# Date: Aug 2026, updated Sep 2026

library(here)
setwd(here::here())

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_sir_parts2_3_single_n.R <nsuccesses>")
}

n <- as.numeric(args[1])
max_budget <- n  # Use n as the computational budget for Parts 3 and 4
cat(sprintf("Running Parts 1, 2, 3, and 4 for nsuccesses = %d\n", n))
cat(sprintf("Parts 3 and 4 will use max_budget = %d\n\n", max_budget))

# Track timing and results
results_summary <- data.frame(
  nsuccesses = integer(),
  part = character(),
  status = character(),
  duration_minutes = numeric(),
  stringsAsFactors = FALSE
)

cat(sprintf("========================================\n"))
cat(sprintf("Processing nsuccesses = %d\n", n))
cat(sprintf("========================================\n\n"))

# # --- Part 1: SIR Maps (Inverse Mapping) ---
# cat(sprintf("[%s] Starting Part 1 (SIR Maps) with n=%d...\n", Sys.time(), n))
# start_time <- Sys.time()
#
# ret_code_1 <- system2(
#   command = "Rscript",
#   args = c(
#     here::here("R", "run_sir_part1_inverse_mapping.R"),
#     as.character(n)
#   ),
#   stdout = here::here("logfiles", sprintf("part1_n%d.Rout", n)),
#   stderr = here::here("logfiles", sprintf("part1_n%d.Rout", n))
# )
#
# end_time <- Sys.time()
# duration_1 <- as.numeric(difftime(end_time, start_time, units = "mins"))
#
# status_1 <- if (ret_code_1 == 0) "SUCCESS" else "FAILED"
# cat(sprintf("[%s] Part 1 completed with status: %s (%.2f minutes)\n",
#             Sys.time(), status_1, duration_1))
#
# results_summary <- rbind(results_summary, data.frame(
#   nsuccesses = n,
#   part = "Part1_SIRMaps",
#   status = status_1,
#   duration_minutes = duration_1,
#   stringsAsFactors = FALSE
# ))
#
# # Save unique timing file for Part 1
# timing_file_part1 <- here::here("data", sprintf("timing_part1_SIRMaps_n%d.csv", n))
# write.csv(data.frame(
#   nsuccesses = n,
#   part = "Part1_SIRMaps",
#   status = status_1,
#   duration_minutes = duration_1,
#   start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
#   end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
#   stringsAsFactors = FALSE
# ), file = timing_file_part1, row.names = FALSE)
# cat(sprintf("Part 1 timing saved to: %s\n\n", timing_file_part1))
#
# # TEMPORARY: Stop after Part 1 (remove this later)
# # stop("Stopping after Part 1 as requested - Parts 2, 3, 4 data already exist")

# # --- Part 2: Basic LHS ---
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

# --- Part 3: OSFD Basic ---
cat(sprintf("[%s] Starting Part 3 (OSFD Basic) with max_budget=%d...\n", Sys.time(), max_budget))
start_time <- Sys.time()

ret_code_3 <- system2(
  command = "Rscript",
  args = c(
    here::here("R", "run_sir_part3_osfd_basic.R"),
    as.character(n),
    as.character(max_budget)  # Pass budget as second argument
  ),
  stdout = here::here("logfiles", sprintf("part3_budget%d.Rout", max_budget)),
  stderr = here::here("logfiles", sprintf("part3_budget%d.Rout", max_budget))
)

end_time <- Sys.time()
duration_3 <- as.numeric(difftime(end_time, start_time, units = "mins"))

status_3 <- if (ret_code_3 == 0) "SUCCESS" else "FAILED"
cat(sprintf("[%s] Part 3 completed with status: %s (%.2f minutes)\n",
            Sys.time(), status_3, duration_3))

results_summary <- rbind(results_summary, data.frame(
  nsuccesses = n,
  part = "Part3_OSFD_Basic",
  status = status_3,
  duration_minutes = duration_3,
  stringsAsFactors = FALSE
))

# Save unique timing file for Part 3
timing_file_part3 <- here::here("data", sprintf("timing_part3_OSFDBasic_n%d.csv", n))
write.csv(data.frame(
  nsuccesses = n,
  part = "Part3_OSFD_Basic",
  status = status_3,
  duration_minutes = duration_3,
  start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
  end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
), file = timing_file_part3, row.names = FALSE)
cat(sprintf("Part 3 timing saved to: %s\n\n", timing_file_part3))

# --- Part 4: Wang Original OSFD ---
cat(sprintf("[%s] Starting Part 4 (Wang Original OSFD) with max_budget=%d...\n", Sys.time(), max_budget))
start_time <- Sys.time()

ret_code_4 <- system2(
  command = "Rscript",
  args = c(
    here::here("R", "run_sir_WangOrig_osfd_basic.R"),
    as.character(n),
    as.character(max_budget)  # Pass budget as second argument
  ),
  stdout = here::here("logfiles", sprintf("part4_budget%d.Rout", max_budget)),
  stderr = here::here("logfiles", sprintf("part4_budget%d.Rout", max_budget))
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

# Save unique timing file for Part 4
timing_file_part4 <- here::here("data", sprintf("timing_part4_WangOrig_n%d.csv", n))
write.csv(data.frame(
  nsuccesses = n,
  part = "Part4_WangOrig",
  status = status_4,
  duration_minutes = duration_4,
  start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
  end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
), file = timing_file_part4, row.names = FALSE)
cat(sprintf("Part 4 timing saved to: %s\n\n", timing_file_part4))

# Print summary for this n
cat(sprintf("\nCompleted n=%d: Part3=%s (%.1f min), Part4=%s (%.1f min)\n",
            n, status_3, duration_3, status_4, duration_4))

cat("\n========================================\n")
cat(sprintf("Task complete for n=%d\n", n))
cat("========================================\n\n")

# Print results
print(results_summary)

# Save individual results
save(results_summary, file = here::here("data", sprintf("sir_parts3_4_n%d_timing.RData", n)))
write.csv(results_summary, file = here::here("data", sprintf("sir_parts3_4_n%d_timing.csv", n)), row.names = FALSE)

cat(sprintf("\nResults saved to:\n"))
cat(sprintf("  - data/sir_parts3_4_n%d_timing.RData\n", n))
cat(sprintf("  - data/sir_parts3_4_n%d_timing.csv\n", n))

# Exit with appropriate code (success only if both parts succeeded)
if (status_3 == "SUCCESS" && status_4 == "SUCCESS") {
  quit(status = 0)
} else {
  quit(status = 1)
}
