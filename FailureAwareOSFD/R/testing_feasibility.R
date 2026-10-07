library(lhs)
library(ggplot2)
library(parallel)

source(here::here("R", "sir_piv_pit_feasible.R"))

# Bounds from run_sir_space_filling.R
pit_bounds <- c(3, 50)
piv_bounds <- c(0.005, 1)
s0 <- 0.999  # Fixed s0 value (within s0_bounds = c(0.95, 0.999))
i0 <- 0.0001 # Initial infected proportion

# Sample using LHS
n_samples <- 10000
lhs_samples <- randomLHS(n_samples, 2)

# Scale to bounds
PIV_samples <- piv_bounds[1] + lhs_samples[, 1] * diff(piv_bounds)
PIT_samples <- pit_bounds[1] + lhs_samples[, 2] * diff(pit_bounds)

# Test feasibility for each point (parallel)
cat("Testing feasibility for", n_samples, "points using 90 cores...\n")
feasible <- unlist(mclapply(1:n_samples, function(i) {
  sir_piv_pit_feasible(PIV_samples[i], PIT_samples[i], s0, i0)
}, mc.cores = 90))

# Create data frame for plotting
plot_data <- data.frame(
  PIV = PIV_samples,
  PIT = PIT_samples,
  feasible = feasible
)

# Summary statistics
cat(sprintf("Feasible points: %d (%.1f%%)\n",
            sum(feasible),
            100 * sum(feasible) / n_samples))

# Create the plot
p <- ggplot(plot_data, aes(x = PIV, y = PIT, color = feasible)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(
    values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
    labels = c("TRUE" = "Feasible", "FALSE" = "Infeasible")
  ) +
  labs(
    title = "SIR Feasibility in PIV-PIT Space",
    subtitle = sprintf("s0 = %.4f, i0 = %.4f", s0, i0),
    x = "Peak Incidence Value (PIV)",
    y = "Peak Incidence Time (PIT)",
    color = "Status"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "top",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Save the plot
ggsave(
  filename = here::here("viz", "sir_feasibility_pivpit.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Plot saved to viz/sir_feasibility_pivpit.png\n")

# # Also create a version with facets for different s0 values
# cat("\nCreating multi-panel plot for different s0 values...\n")
# s0_values <- c(0.95, 0.97, 0.99, 0.999)
# multi_plot_data <- data.frame()

# for(s0_val in s0_values) {
#   i0_val <- 1 - s0_val

#   # Test feasibility
#   feasible_s0 <- sapply(1:n_samples, function(i) {
#     sir_piv_pit_feasible(PIV_samples[i], PIT_samples[i], s0_val, i0_val)
#   })

#   # Append to data frame
#   temp_df <- data.frame(
#     PIV = PIV_samples,
#     PIT = PIT_samples,
#     feasible = feasible_s0,
#     s0 = sprintf("s0 = %.3f", s0_val)
#   )

#   multi_plot_data <- rbind(multi_plot_data, temp_df)

#   cat(sprintf("  s0 = %.3f: %d feasible (%.1f%%)\n",
#               s0_val,
#               sum(feasible_s0),
#               100 * sum(feasible_s0) / n_samples))
# }

# # Create multi-panel plot
# p_multi <- ggplot(multi_plot_data, aes(x = PIV, y = PIT, color = feasible)) +
#   geom_point(alpha = 0.4, size = 0.8) +
#   facet_wrap(~ s0, ncol = 2) +
#   scale_color_manual(
#     values = c("TRUE" = "#2E7D32", "FALSE" = "#C62828"),
#     labels = c("TRUE" = "Feasible", "FALSE" = "Infeasible")
#   ) +
#   labs(
#     title = "SIR Feasibility in PIV-PIT Space Across s0 Values",
#     x = "Peak Incidence Value (PIV)",
#     y = "Peak Incidence Time (PIT)",
#     color = "Status"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(face = "bold", size = 14),
#     legend.position = "top",
#     strip.text = element_text(face = "bold")
#   )

# # Save multi-panel plot
# ggsave(
#   filename = here::here("viz", "sir_feasibility_pivpit_multi_s0.png"),
#   plot = p_multi,
#   width = 10,
#   height = 8,
#   dpi = 300
# )

# cat("Multi-panel plot saved to viz/sir_feasibility_pivpit_multi_s0.png\n")