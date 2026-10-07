# Plots for output space filling designs by regime.
## Author: AC Murph
## Date February 2026

# Script to plot the SIR curve space-filling experiment.
# Author: AC Murph
library(ggplot2)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(latex2exp)
library(reshape2)


####################
### Prepare sMOA ###
####################
k <- 5
h <- 4
closest <- 4422

# Load pre-computed synthetic embeddings (same pattern as perform_smoa_using_new_synthetic.R)
embeddings_file <- here::here('data', 'synthetic_embeddings.RData')
if (!file.exists(embeddings_file)) {
  stop("Synthetic embeddings not found. Run R/prepare_synthetic_embeddings.R first.")
}
load(embeddings_file)
print(paste("Loaded embeddings from", embeddings_file))
print(paste("ISFD dimensions: X =", paste(dim(embed_mat_X_ISFD), collapse=" x "),
            ", y =", paste(dim(embed_mat_y_ISFD), collapse=" x ")))
print(paste("OSFD dimensions: X =", paste(dim(embed_mat_X_OSFD), collapse=" x "),
            ", y =", paste(dim(embed_mat_y_OSFD), collapse=" x ")))

# Combine ISFD and OSFD data with method and regime labels
embed_mat_X_ISFD_df <- as.data.frame(embed_mat_X_ISFD)
colnames(embed_mat_X_ISFD_df) <- as.character(1:k)
embed_mat_X_ISFD_df$method <- "ISFD"
embed_mat_X_ISFD_df$regime <- regime_ISFD

embed_mat_X_OSFD_df <- as.data.frame(embed_mat_X_OSFD)
colnames(embed_mat_X_OSFD_df) <- as.character(1:k)
embed_mat_X_OSFD_df$method <- "OSFD"
embed_mat_X_OSFD_df$regime <- regime_OSFD

# Combine both
embed_mat_X <- rbind(embed_mat_X_ISFD_df, embed_mat_X_OSFD_df)
embed_mat_X$sim_num <- 1:nrow(embed_mat_X)

# Reshape for plotting
input_space <- reshape2::melt(embed_mat_X, id.vars = c("method", "regime", "sim_num"))
n_sims_per_group <- 200  

sims_keep <- input_space %>%
  distinct(method, regime, sim_num) %>%
  group_by(method, regime) %>%
  slice_sample(
    n = n_sims_per_group,   # guard if a group has < n_sims_per_group sims
    replace = FALSE
    # alternatively:
    # prop = prop_sims_per_group
  ) %>%
  ungroup()

# 2) keep *all* rows matching those sampled sims
input_space_sub <- input_space %>%
  semi_join(sims_keep, by = c("method", "regime", "sim_num"))
# --- plot (free y across facets) ---
p = ggplot(input_space_sub, aes(x = variable, y = value, color = method, group = sim_num)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ regime, nrow = 2, ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(x = "Snippet Point Number", y = "Value", color = "method") +
  theme(
    plot.title = element_text(size = 40),     # not bold
    axis.title = element_text(size = 34),     # not bold
    axis.text  = element_text(size = 26),

    strip.text = element_text(size = 34),     # facet labels bigger
    # strip.background = element_rect(linewidth = 1.2),  # optional

    legend.title = element_text(size = 34),   # legend title bigger
    legend.text  = element_text(size = 28),   # legend entries bigger
    legend.key.size = unit(1.6, "cm"),        # bigger legend keys
    legend.spacing.y = unit(0.6, "cm")        # more breathing room (optional)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

scale = 0.7
ggsave(
  filename = paste0(here::here("viz", "sir_by_regime_experiment.png")),
  plot = p,
  width = 40*scale,      # adjust width as needed
  height = 20*scale,     # adjust height as needed
  dpi = 300       # high-quality resolution
)

