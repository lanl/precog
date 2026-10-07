# Create work manifest for parallel processing
## This script creates a CSV mapping of all (disease, time_series_index) pairs
## to enable even distribution of work across SLURM nodes
## Author: LJ Beesley + AC Murph
## Date: August 2026

library(here)

setwd(here::here())

data_path <- here::here("data", "raw_data")
output_path <- here::here("data")

# Get the same file list as in the main script
FILES_ALL <- gsub('.RDS', '', list.files(data_path))
FILES_ALL <- FILES_ALL[!grepl('cimmid', FILES_ALL)]
FILES_ALL <- FILES_ALL[!grepl('ginkgo', FILES_ALL)]
FILES_ALL <- FILES_ALL[grepl('Chikungunya_deSouza', FILES_ALL) |
                       grepl('Dengue_opendengue', FILES_ALL) |
                       grepl('Influenza_usflunet', FILES_ALL) |
                       grepl('Influenza_ushhs', FILES_ALL) |
                       grepl('Mpox_who', FILES_ALL) |
                       grepl('COVID', FILES_ALL)]

cat(sprintf("Found %d disease datasets\n", length(FILES_ALL)))

# Create manifest of all (disease, time_series_index) pairs
manifest <- data.frame()

for (disease in FILES_ALL) {
  cat(sprintf("Processing %s...\n", disease))

  # Load the data
  list_of_lists <- readRDS(file.path(data_path, paste0(disease, ".RDS")))
  n_ts <- length(list_of_lists)

  cat(sprintf("  Found %d time series\n", n_ts))

  # Add to manifest
  disease_manifest <- data.frame(
    disease = disease,
    ts_index = 1:n_ts
  )

  manifest <- rbind(manifest, disease_manifest)
}

# Add a sequential work_id
manifest$work_id <- 1:nrow(manifest)

cat(sprintf("\nTotal work items: %d\n", nrow(manifest)))
cat(sprintf("Items per node (if using 6 nodes): %.1f\n", nrow(manifest) / 6))

# Save manifest
manifest_file <- file.path(output_path, "work_manifest.csv")
write.csv(manifest, manifest_file, row.names = FALSE, quote = FALSE)
cat(sprintf("\nManifest saved to: %s\n", manifest_file))

# Print summary by disease
cat("\nBreakdown by disease:\n")
summary_table <- table(manifest$disease)
print(summary_table)
