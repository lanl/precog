# Quick check of the feasible X_isfd samples
library(ggplot2)
library(patchwork)

# Load the workspace from run_sir_space_filling.R
# Assumes X_isfd, bounds, etc. are in environment

# Scale X_isfd back to original bounds
PIV_scaled <- piv_bounds[1] + X_isfd[, 1] * diff(piv_bounds)
PIT_scaled <- pit_bounds[1] + X_isfd[, 2] * diff(pit_bounds)
s0_scaled <- s0_bounds[1] + X_isfd[, 3] * diff(s0_bounds)

# Create data frame for plotting
plot_df <- data.frame(
  PIV = PIV_scaled,
  PIT = PIT_scaled,
  s0 = s0_scaled
)

cat(sprintf("X_isfd contains %d feasible samples\n", nrow(X_isfd)))
cat(sprintf("PIV range: [%.4f, %.4f]\n", min(PIV_scaled), max(PIV_scaled)))
cat(sprintf("PIT range: [%.1f, %.1f]\n", min(PIT_scaled), max(PIT_scaled)))
cat(sprintf("s0 range:  [%.4f, %.4f]\n", min(s0_scaled), max(s0_scaled)))

# Main plot: PIV vs PIT
p1 <- ggplot(plot_df, aes(x = PIV, y = PIT)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#1976D2") +
  geom_density_2d(color = "#D32F2F", alpha = 0.7) +
  labs(
    title = "Feasible PIV-PIT Space",
    subtitle = sprintf("n = %d feasible samples from LHS filtering", nrow(plot_df)),
    x = "Peak Incidence Value (PIV)",
    y = "Peak Incidence Time (PIT)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

# PIV vs s0
p2 <- ggplot(plot_df, aes(x = PIV, y = s0)) +
  geom_point(alpha = 0.5, size = 1, color = "#388E3C") +
  labs(
    title = "PIV vs s0",
    x = "Peak Incidence Value (PIV)",
    y = "Initial Susceptible (s0)"
  ) +
  theme_minimal()

# PIT vs s0
p3 <- ggplot(plot_df, aes(x = PIT, y = s0)) +
  geom_point(alpha = 0.5, size = 1, color = "#7B1FA2") +
  labs(
    title = "PIT vs s0",
    x = "Peak Incidence Time (PIT)",
    y = "Initial Susceptible (s0)"
  ) +
  theme_minimal()

# Histograms
p4 <- ggplot(plot_df, aes(x = PIV)) +
  geom_histogram(bins = 50, fill = "#1976D2", alpha = 0.7) +
  labs(title = "PIV Distribution", x = "PIV", y = "Count") +
  theme_minimal()

p5 <- ggplot(plot_df, aes(x = PIT)) +
  geom_histogram(bins = 50, fill = "#D32F2F", alpha = 0.7) +
  labs(title = "PIT Distribution", x = "PIT", y = "Count") +
  theme_minimal()

p6 <- ggplot(plot_df, aes(x = s0)) +
  geom_histogram(bins = 50, fill = "#388E3C", alpha = 0.7) +
  labs(title = "s0 Distribution", x = "s0", y = "Count") +
  theme_minimal()

# Combine plots
p_combined <- (p1 | (p2 / p3)) / (p4 | p5 | p6) +
  plot_annotation(
    title = "Feasible Space-Filling Sample Visualization",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

# Save
ggsave(
  filename = here::here("viz", "feasible_samples_check.png"),
  plot = p_combined,
  width = 14,
  height = 10,
  dpi = 300
)

cat("\nPlot saved to viz/feasible_samples_check.png\n")

# Show main plot
print(p1)
