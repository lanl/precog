#rm(list=ls())

#############################
### Plots for publication ###
#############################

# basic libraries
library(rstudioapi)
library(geomtextpath)
library(dplyr)  
library(lubridate)
library(data.table)
library(effectsize)
library(tidyr)
library(rcompanion)
library(boot)

# plotting libraries
library(ggplot2) 
library(GGally) 
library(maps) 
library(viridis)
library(cowplot)
library(gridExtra)
library(ggExtra)
library(latex2exp)
library(RColorBrewer)

dir = "./precog/reporting-delay/"
results_dir = paste0(dir, "results/")
output_dir = paste0(dir, "plots/")

##############################
### Read in data and clean ###
##############################

### GISAID data for data section figure ###
merged_df <- fread(paste0(dir, "reporting_delay_data.csv")) %>%
  select(-V1) %>% 
  rename(Date = "collection_date",
         Location = "Admin0") %>%
  mutate(Year = year(Date)) %>%
  filter(Location %in% c("Global", "United States", "Denmark", "Brazil")) %>%
  mutate(Location = factor(Location, levels = c("Brazil", "Denmark", "United States", "Global")))

### metric results for all locations ###
all_locations <- fread(paste0(results_dir, "all_metric_results.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2020-11-01")) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant") %>%
  select(-c(V1, error_code)) %>%
  arrange(Date) %>%
  mutate(r_cat = cut(r, 
                         breaks = seq(0, 100, 10), 
                         labels = c("(0, 10]", "(10, 20]", "(20, 30]", "(30, 40]",
                                    "(40, 50]", "(50, 60]", "(60, 70]", "(70, 80]",
                                    "(80, 90]", "(90, 100)"), 
                         include.lowest = T),
         n_int = cut(n,
                     breaks = c(0, 100, 1000, 10000, Inf),
                     labels = c("n: 2-100", "n: 101-1,000", "n: 1,001-10,000", 
                                "n > 10,000"),
                     include.lowest = T)
  )

### simulations under the null ###
load(file = paste0(results_dir, "all_sim_results.RData")) # this is metric_comb

### results for emerging variant time series ###
emerge_join <- fread(paste0(results_dir, "emerging_variant_ts.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  select(-c(V1)) %>%
  left_join(metric_comb)
  
####################
### Data section ###
####################

### reporting delay histogram by location and year ###
merged_df %>%
  filter(delay_days <= 10*7,
         Year < 2023
         ) %>%
  mutate(delay_weeks = delay_days/7) %>%
  ggplot(aes(x = delay_weeks, weight = counts, 
             fill = Location)) +
  geom_density(color = "black", alpha = 0.5, show.legend = T) +
  scale_fill_manual(values = brewer.pal(4, name = "Dark2")) +
  labs(x = "Reporting Delay (Weeks)",
       y = "Sequence Density") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        panel.spacing.x = unit(2, "lines"),
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        plot.margin = margin(5.5, 25, 5.5, 5.5)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(0, 10, 2)) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ Year, ncol = 3, scales = "free_y")
ggsave(filename = "compare_delays.png", path = output_dir,
       width = 12, height = 4, unit = "in", limitsize = F)

### percent reported within 1 and 2 weeks ###

# disaggregated by year
merged_df %>%
  filter(Year < 2023) %>%
  group_by(Location, Year) %>%
  summarize(n_total = sum(counts, na.rm = T),
            n_1wk = sum(counts[delay_days <= 7], na.rm = T),
            pct_1wk = 100 * (n_1wk/n_total),
            n_2wk = sum(counts[delay_days <= 14], na.rm = T),
            pct_2wk = 100 * (n_2wk/n_total),
            .groups = "drop") %>%
  select(-n_1wk, -n_2wk)

# entire time frame
merged_df %>%
  filter(Year < 2023) %>%
  group_by(Location) %>%
  summarize(n_total = sum(counts, na.rm = T),
            n_1wk = sum(counts[delay_days <= 7], na.rm = T),
            pct_1wk = 100 * (n_1wk/n_total),
            n_2wk = sum(counts[delay_days <= 14], na.rm = T),
            pct_2wk = 100 * (n_2wk/n_total),
            .groups = "drop") %>%
  select(-n_1wk, -n_2wk)

######################
### n-r categories ###
######################

# categories align with w_vs_reporting_rate.png

### all countries and delay periods ###
n_r_cats <- metric_comb %>% 
  filter(n > 1, r < 100) %>%
  mutate(n_cat = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("(1, 50]", "(50, 500]", 
                                "(500, 5k]", "> 5k"), 
                     include.lowest = T)) %>%
  count(n_cat, r_cat) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = r_cat, y = n_cat, fill = prop*100)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1),
                color = prop > 0.15), size = 5) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "Sample Size (n)") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Percent of Observations",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")


# statistical significance
n_r_sig <- metric_comb %>% 
  filter(n > 1, r < 100,
         Metric == "w") %>%
  mutate(n_cat = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("(1, 50]", "(50, 500]", 
                                "(500, 5k]", "> 5k"), 
                     include.lowest = T)) %>%
  count(n_cat, r_cat, pvalue_cat) %>%
  pivot_wider(names_from = pvalue_cat, values_from = n) %>%
  mutate(sig_prop = `< 0.05` / (`< 0.05` + `> 0.05`),
         sig_prop = case_when(is.na(`< 0.05`) ~ 0,
                              is.na(`> 0.05`) ~ 1,
                              .default = sig_prop)) %>% 
  ggplot(aes(x = r_cat, y = n_cat, fill = sig_prop)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(sig_prop, accuracy = 0.1),
                color = sig_prop > 0.2), size = 5) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Significance \nProportion",
                       na.value = "white", option = "magma",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")


plot_grid(n_r_cats, n_r_sig, ncol = 2, align = "h", axis = "l")
ggsave(filename = "n_r_cats_all_countries.png", path = output_dir,
       width = 16, height = 6, unit = "in")

### separate delay periods ###
all_countries_cat <- metric_comb %>% 
  filter(n > 1,
         Location != "Global",
         Metric == "w") %>%
  mutate(r_cat = cut(r,
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c("(0, 20]", "(20, 40]", "(40, 60]",
                               "(60, 80]", "(80, 100]"),
                     include.lowest = T),
         n_cat = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("(1, 50]", "(50, 500]", 
                                "(500, 5k]", "> 5k"), 
                     include.lowest = T)) %>%
  group_by(delay_days) %>%
  count(n_cat, r_cat) %>%
  mutate(prop = n/ sum(n),
         Location = "All Countries") %>%
  select(Location, delay_days, n_cat, r_cat, n, prop)


by_location_cat <- metric_comb %>% 
  filter(n > 0, #omega < 100, 
         Metric == "w",
         Location %in% c("Brazil", "Denmark", "United States")) %>%
  mutate(r_cat = cut(r,
                     breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c("(0, 20]", "(20, 40]", "(40, 60]",
                                "(60, 80]", "(80, 100]"),
                     include.lowest = T),
         n_cat = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("(1, 50]", "(50, 500]", 
                                "(500, 5k]", "> 5k"), 
                     include.lowest = T)) %>%
  group_by(Location, delay_days) %>%
  count(n_cat, r_cat) %>%
  mutate(prop = n/ sum(n))
  
# combine for plot
rbind(all_countries_cat, by_location_cat) %>%
  ggplot(aes(x = r_cat, y = n_cat, fill = prop*100)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1),
                color = prop > 0.3), size = 3) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "Sample Size (n)") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Percent of Observations",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_grid(Location ~ delay_days)
ggsave(filename = "n_r_cats.png", path = output_dir,
       width = 10, height = 10, unit = "in")

##########################  
### Cohen's w vs Norms ###
##########################

## when considering observations with at least 2 near-real-time samples
#filter_n <- all_locations %>% filter(n > 1)
filter_n <- all_locations

cor(filter_n$w, filter_n$L1,
    use = 'pairwise.complete.obs') %>% round(3)
cor(filter_n$w, filter_n$L2,
    use = 'pairwise.complete.obs') %>% round(3)
cor(filter_n$w, filter_n$Linf,
    use = 'pairwise.complete.obs') %>% round(3)

cor(filter_n$w, filter_n$L1,
    use = 'pairwise.complete.obs',
    method = 'spearman') %>% round(3)
cor(filter_n$w, filter_n$L2,
    use = 'pairwise.complete.obs',
    method = 'spearman') %>% round(3)
cor(filter_n$w, filter_n$Linf,
    use = 'pairwise.complete.obs',
    method = 'spearman') %>% round(3)

### boxplot ###
p <- all_locations %>%
  filter(n > 1,
         !(is.na(w))) %>%
  select(w, L1, L2, Linf, r) %>%
  pivot_longer(cols = c(L1, L2, Linf), names_to = "Norm", values_to = "Norm_Value") %>%
  mutate(norm_int = cut(Norm_Value,
                        breaks = seq(0, 2, 0.25),
                        labels = c("(0, 0.25]", "(0.25, 0.5]", "(0.5, 0.75]", "(0.75, 1]",
                                   "(1, 1.25]", "(1.25, 1.5]", "(1.5, 1.75]", "(1.75, 2]"),
                        include.lowest = T)
  ) %>%
  ggplot(aes(x = norm_int, y = w, color = r)) +
  geom_point(position = position_jitterdodge(), alpha = 0.2, size = 1, show.legend = T) +
  geom_boxplot(alpha = 0) +
  theme_classic() +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "Reporting Rate (%)") +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = 20, barheight = 1)) +
  labs(x = "Norm", 
       y = "Cohen\'s w") +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  coord_cartesian(ylim=c(0,10)) +
  facet_wrap(~ Norm, ncol = 1)
ggsave(filename = "w_vs_Norms_box.png", p,
       path = output_dir, width = 7, height = 8, unit = "in")


#######################
### Metrics vs n, r ###
#######################

### Cohen's w ###
metric_comb %>%
  filter(Metric == "w",
         r < 100,
         n > 1) %>%
  mutate(
    n_int = cut(n,
                breaks = c(0, 50, 500, 5000, Inf),
                labels = c("n: 2-50", "n: 51-500", "n: 501-5,000", "n > 5,000"),
                include.lowest = T)
         ) %>%
  ggplot() +
  geom_jitter(aes(x = factor(r_cat), y = True_Value, color = pvalue_cat),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(r_cat), y = True_Value),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "Cohen's w") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Set1"),
                     name = "p-value") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ n_int, nrow = 1) +
  coord_cartesian(ylim = c(0, 1))
ggsave(filename = "w_vs_reporting_rate.png", path = output_dir,
       width = 15, height = 5, unit = "in")

denom <- metric_comb %>%
  filter(Metric == "V",
         omega < 100,
         n > 1,
         #n <= 5000,
         True_Value >= 1,
         )
num <- denom %>% filter(pvalue < 0.05)

nrow(num)/nrow(denom)  

### L1 ###
all_locations %>%
  filter(n > 1) %>%
  mutate(
    n_int = cut(n,
                breaks = c(0, 50, 500, 5000, Inf),
                labels = c("n: 2-50", "n: 51-500", "n: 501-5,000", "n > 5,000"),
                include.lowest = T)
  ) %>%
  ggplot() +
  geom_jitter(aes(x = factor(r_cat), y = L1, color = w),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(r_cat), y = L1),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "L1 Norm") +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "w",
                      limits = c(0, 1),
                      na.value = viridis(n = 1)
                      ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 2)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ n_int, nrow = 1)
ggsave(filename = "L1_vs_reporting_rate.png", path = output_dir,
       width = 15, height = 5, unit = "in")

### Fei ###
all_fei <- all_locations %>%
  filter(n > 1) %>%
  mutate(
    n_int = cut(n,
                breaks = c(0, 50, 500, 5000, Inf),
                labels = c("n: 2-50", "n: 51-500", "n: 501-5,000", "n > 5,000"),
                include.lowest = T)
  ) %>%
  ggplot() +
  geom_jitter(aes(x = factor(r_cat), y = Fei, color = w),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(r_cat), y = Fei),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "Fei") +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "w",
                      limits = c(0, 1),
                      na.value = viridis(n = 1)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ n_int, nrow = 1)

fei_0.1 <- all_locations %>%
  filter(n > 1) %>%
  mutate(
    n_int = cut(n,
                breaks = c(0, 50, 500, 5000, Inf),
                labels = c("n: 2-50", "n: 51-500", "n: 501-5,000", "n > 5,000"),
                include.lowest = T)
  ) %>%
  ggplot() +
  geom_jitter(aes(x = factor(r_cat), y = Fei, color = w),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(r_cat), y = Fei),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "Fei") +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "w",
                      limits = c(0, 1),
                      na.value = viridis(n = 1)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ n_int, nrow = 1) +
  coord_cartesian(ylim = c(0, 0.1))

fei_0.015 <- all_locations %>%
  filter(n > 1) %>%
  mutate(
    n_int = cut(n,
                breaks = c(0, 50, 500, 5000, Inf),
                labels = c("n: 2-50", "n: 51-500", "n: 501-5,000", "n > 5,000"),
                include.lowest = T)
  ) %>%
  ggplot() +
  geom_jitter(aes(x = factor(r_cat), y = Fei, color = w),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(r_cat), y = Fei),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "Fei") +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "w",
                      limits = c(0, 1),
                      na.value = viridis(n = 1)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ n_int, nrow = 1) +
  coord_cartesian(ylim = c(0, 0.015))

plot_grid(all_fei, fei_0.1, fei_0.015, ncol = 1, align = "v", axis = "l")
ggsave(filename = "Fei_vs_reporting_rate.png", path = output_dir,
       width = 15, height = 15, unit = "in")


####################################
### Simulations - p-values & q95 ###
####################################

### q95 vs n ###
q95_vs_n <- metric_comb %>%
  filter(n > 1, r < 100) %>%
  mutate(q95 = ifelse(q95 < 1e-4, 1e-4, q95)) %>%
  group_by(n_cat, Metric) %>%
  summarize(med = median(q95, na.rm = T),
            quant_05 = quantile(q95, probs = 0.25, na.rm = T),
            quant_95 = quantile(q95, probs = 0.75, na.rm = T)) %>%
  ggplot() +
  geom_ribbon(aes(x = as.numeric(n_cat), ymin = quant_05, ymax = quant_95, fill = Metric),
              alpha = 0.15) +
  geom_line(aes(x = as.numeric(n_cat), y = med, color = Metric), lwd = 2) +
  theme_classic() +
  labs(x = "Sample Size (n)", 
       y = "95th Percentile",
       title = "") +
  scale_color_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:length(levels(metric_comb$n_cat)),
                     labels = c("(1, 10]", "(10, 100]", "(100, 1k]", 
                                "(1k, 10k]", "> 10k"),
                     expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")


### q95 vs r ###
q95_vs_r <- metric_comb %>%
  filter(n > 1, r < 100) %>%
  mutate(q95 = ifelse(q95 < 1e-4, 1e-4, q95)) %>%
  group_by(r_cat, Metric) %>%
  summarize(med = median(q95, na.rm = T),
            quant_05 = quantile(q95, probs = 0.25, na.rm = T),
            quant_95 = quantile(q95, probs = 0.75, na.rm = T)) %>%
  ggplot() +
  geom_ribbon(aes(x = as.numeric(r_cat), ymin = quant_05, ymax = quant_95, fill = Metric),
              alpha = 0.15) +
  geom_line(aes(x = as.numeric(r_cat), y = med, color = Metric), lwd = 2) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "95th Percentile",
       title = "") +
  scale_color_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:length(levels(metric_comb$r_cat)),
                     labels = levels(metric_comb$r_cat),
                     expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.spacing = unit(2, "cm"))


plot_grid(q95_vs_n, q95_vs_r, ncol = 2, align = "hv", axis = "lr")
ggsave(filename = "q95.png", path = output_dir,
       width = 12, height = 6, unit = "in")

########################################
### Simulations - Metric Time Series ###
########################################

# find when variants emerged by location:
#   Alpha, Delta, Omicron BA.1, Omicron BA.2, Omicron BA.4/5

find_emergence <- function(loc, variant, min_prop = 0.05){
  if(variant != "Omicron BA.4"){
    variant_subset <- emerge_join %>%
      filter(Location == loc,
             pango == variant,
             p_val >= min_prop)
  } else {
    variant_subset <- emerge_join %>% 
      filter(Location == loc, 
             pango %in% c("Omicron BA.4", "Omicron BA.5"), 
             Metric == "w") %>% 
      group_by(Date) %>% 
      summarize(sum_prop = sum(p_val)) %>% 
      filter(sum_prop >= min_prop)
  }
  return(min(variant_subset$Date))
}  

global_alpha <- find_emergence("Global", "Alpha")
brazil_alpha <- find_emergence("Brazil", "Alpha")
denmark_alpha <- find_emergence("Denmark", "Alpha")
us_alpha <- find_emergence("United States", "Alpha")

brazil_alpha - global_alpha
denmark_alpha - global_alpha
us_alpha - global_alpha

global_delta <- find_emergence("Global", "Delta")
brazil_delta <- find_emergence("Brazil", "Delta")
denmark_delta <- find_emergence("Denmark", "Delta")
us_delta <- find_emergence("United States", "Delta")

brazil_delta - global_delta
denmark_delta - global_delta
us_delta - global_delta

global_BA.1 <- find_emergence("Global", "Omicron BA.1")
brazil_BA.1 <- find_emergence("Brazil", "Omicron BA.1")
denmark_BA.1 <- find_emergence("Denmark", "Omicron BA.1")
us_BA.1 <- find_emergence("United States", "Omicron BA.1")

brazil_BA.1 - global_BA.1
denmark_BA.1 - global_BA.1
us_BA.1 - global_BA.1

global_BA.2 <- find_emergence("Global", "Omicron BA.2")
brazil_BA.2 <- find_emergence("Brazil", "Omicron BA.2")
denmark_BA.2 <- find_emergence("Denmark", "Omicron BA.2")
us_BA.2 <- find_emergence("United States", "Omicron BA.2")

brazil_BA.2 - global_BA.2
denmark_BA.2 - global_BA.2
us_BA.2 - global_BA.2

global_BA.4 <- find_emergence("Global", "Omicron BA.4")
brazil_BA.4 <- find_emergence("Brazil", "Omicron BA.4")
denmark_BA.4 <- find_emergence("Denmark", "Omicron BA.4")
us_BA.4 <- find_emergence("United States", "Omicron BA.4")

brazil_BA.4 - global_BA.4
denmark_BA.4 - global_BA.4
us_BA.4 - global_BA.4

variant_dates <- c(global_alpha, global_delta, global_BA.1, global_BA.2, global_BA.4)

# true metric time series, alpha by significance
compare_w <- metric_comb %>%
  filter(Metric == "w",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1
         ) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection")),
         True_Value = ifelse(True_Value > 1, 1, True_Value)
         ) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = True_Value)) +
  # geom_vline(
  #   xintercept = variant_dates,
  #   color = "gray40",
  #   linetype = "dashed",
  #   linewidth = 1.5,
  #   alpha = 0.8
  # ) +
  geom_line(aes(color = Location),
            lwd = 1, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = Location, shape = Location),
            #color = "black", 
            size = 2, show.legend = F) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_fill_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  theme_bw() +
  labs(x = "", 
       y = "Cohen's w",
       title = "") +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  facet_grid(~ delay_days)


# n time series
compare_n <- metric_comb %>%
  filter(Metric == "w",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1
         ) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = n)) +
  # geom_vline(
  #   xintercept = variant_dates,
  #   color = "gray40",
  #   linetype = "dashed",
  #   linewidth = 1.5,
  #   alpha = 0.8
  # ) +
  geom_line(aes(color = Location),
            lwd = 1, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = Location, shape = Location),
             #color= "black", 
             alpha = 0.7, size = 2, show.legend = F) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_fill_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  theme_bw() +
  labs(x = "", 
       y = "Sample Size (n)",
       title = "") +
  scale_y_log10() + 
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  facet_wrap(~ delay_days, ncol = 2)


# r time series
compare_r <- metric_comb %>%
  filter(Metric == "w",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1
         ) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = r)) +
  # geom_vline(
  #   xintercept = variant_dates,
  #   color = "gray40",
  #   linetype = "dashed",
  #   linewidth = 1.5,
  #   alpha = 0.8
  # ) +
  geom_line(aes(color = Location),
            lwd = 1, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = Location, shape = Location),
             #color = "black", 
             alpha = 0.7, size = 2, show.legend = F) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_fill_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  theme_bw() +
  labs(x = "", 
       y = "Reporting Rate (%)",
       title = "") +
  ylim(c(0, 100)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  facet_wrap(~ delay_days, ncol = 2)


# w - q95 time series
w_minus_q95 <- metric_comb %>%
  filter(Metric == "w",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1
         ) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  mutate(dif = True_Value - q95) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = dif)) +
  # geom_vline(
  #   xintercept = variant_dates,
  #   color = "gray40",
  #   linetype = "dashed",
  #   linewidth = 1.5,
  #   alpha = 0.8
  # ) +
  geom_line(aes(color = Location),
            lwd = 1, alpha = 0.5, show.legend = F) +
  geom_point(aes(fill = Location, shape = Location),
             #color = "black", 
             alpha = 0.7, size = 2, show.legend = T) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_fill_manual(values = brewer.pal(4, name = "Dark2")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  theme_bw() +
  labs(x = "Collection Date", 
       y = "w - 95th Percentile",
       title = "") +
  ylim(c(-0.1, 1)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  facet_wrap(~ delay_days, ncol = 2)


plot_grid(compare_w, compare_n, compare_r, w_minus_q95, ncol = 1, align = "v", axis = "l")
ggsave(filename = "Compare_w_Locations.png", path = output_dir,
       width = 16, height = 16, unit = "in")

#######################################
### Emerging Variants - time series ###
#######################################

pango_colors = c(brewer.pal(n = 8, name = "Dark2"), 
                 brewer.pal(n = 11, name = "RdYlBu")[c(2, 4, 8, 11)],
                 brewer.pal(n = 11, name = "PRGn")[c(2, 4, 8, 11)])

first_date = as.Date('2020-11-01')
last_date = as.Date('2022-12-31')

make_emergence_plot <- function(loc, delay, cap_w = FALSE) {
  
  date_seq   <- seq(first_date, last_date, by = "day")
  
  loc_title <- ifelse(loc == "United States", "US", loc)
  
  sub_df <- emerge_join %>%
    filter(Location == loc,
           delay_days == delay,
           Date %in% date_seq) %>%
    arrange(Date)
  
  # Order variants according to the date of their maximum validation proportion
  peak_max <- sub_df %>%
    group_by(pango) %>%
    slice_max(order_by = p_val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(pango, peak_time = Date)
  
  max_props <- sub_df %>%
    group_by(pango) %>%
    summarize(
      max_p_val = max(p_val, na.rm = TRUE),
      max_p_nrt = max(p_nrt, na.rm = TRUE),
      overall_prop_dif = max_p_val - max_p_nrt,
      .groups = "drop"
    ) %>%
    filter(
      max_p_val >= 0.05 | abs(overall_prop_dif) > 0.05) %>%
    left_join(peak_max, by = "pango") %>%
    arrange(peak_time)
  
  sub_df <- sub_df %>% filter(pango %in% max_props$pango)
  
  pango_order <- unique(c("Other/Unknown", max_props$pango))
  
  # Cohen's w
  w_data <- sub_df %>%
    filter(Metric == "w", n > 1) %>%
    mutate(pvalue_cat = factor(pvalue_cat, levels = c("> 0.05", "< 0.05")),
      True_Value = if(cap_w) {pmin(True_Value, 1)} else {True_Value}
    )
  
  w_plot <- ggplot(w_data, aes(x = Date, y = True_Value, color = pvalue_cat)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_color_manual(values = brewer.pal(n = 3, name = "Paired"),
                       name = "p-value",
                       na.translate = FALSE) +
    scale_x_date(limits = c(first_date, last_date), 
                 expand = expansion(mult = 0)) +
    labs(x = "", y = "Cohen's w",
         title = paste0(loc, ": d = ", delay, " Days Post-Collection")) +
    theme_classic() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    )
  
  # Validation proportion
  p_val_data <- sub_df %>%
    mutate(pango = factor(pango, levels = pango_order))
  
  p_val_plot <- ggplot(p_val_data, aes(x = Date, y = p_val, color = pango)) +
    geom_line(linewidth = 1) +
    geom_point(aes(fill = pango), alpha = 0.7, show.legend = FALSE) +
    scale_color_manual(values = pango_colors, name = "", 
                       na.translate = FALSE) +
    scale_x_date(limits = c(first_date, last_date), 
                 expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "", y = "Validation Proportion",
         title = TeX("Validation Proportion (${p}_{val}$)")) +
    theme_classic() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.spacing = unit(5, "cm")
    )
  
  # Proportion difference
  prop_dif_data <- sub_df %>%
    filter(Metric == "w", n > 1) %>%
    mutate(pango = factor(pango, levels = pango_order)) %>%
    group_by(pango) %>%
    complete(Date = date_seq) %>%
    arrange(Date, .by_group = TRUE) %>%
    mutate(line_group = cumsum(is.na(prop_dif))) %>%
    ungroup()
  
  prop_dif_plot <- ggplot(prop_dif_data,
                          aes(x = Date, y = prop_dif, color = pango)) +
    geom_line(aes(group = interaction(pango, line_group)),
              alpha = 0.7, linewidth = 1) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    scale_color_manual(values = pango_colors, name = "",
                       na.translate = FALSE) +
    scale_x_date(limits = c(first_date, last_date),
                 expand = expansion(mult = 0)) +
    labs(x = "Collection Date", y = "Proportion Difference",
         title = TeX("Proportion Difference: ${p}_{val} - \\hat{p}_{nrt}(d)$")) +
    theme_classic() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.spacing = unit(2, "cm")
    )
  
  combined_plot <- plot_grid(w_plot, p_val_plot, prop_dif_plot,
                             ncol = 1, align = "v", axis = "lr")
  
  ggsave(filename = paste(loc_title, delay, "emerge.png", sep = "_"),
         path = output_dir, plot = combined_plot,
         width = 12, height = 12, units = "in")
  
  return(combined_plot)
}

make_emergence_plot(loc = "United States", delay = 7)
make_emergence_plot(loc = "United States", delay = 14)
make_emergence_plot(loc = "Denmark", delay = 7, cap_w = T)
make_emergence_plot(loc = "Denmark", delay = 14)
make_emergence_plot(loc = "Brazil", delay = 7)
make_emergence_plot(loc = "Brazil", delay = 14)


####################################
### sensitivity analysis: vary M ###
####################################

M_1000 <- metric_comb %>%
  filter(n > 1, r < 100) %>%
  select(Metric, q95, n_cat) %>%
  mutate(M = 1000)

load(file = paste0(results_dir, "sim_results_M_500.RData"))
M_500 <- metric_comb %>%
  filter(n > 1, r < 100) %>%
  select(Metric, q95, n_cat) %>%
  mutate(M = 500)

load(file = paste0(results_dir, "sim_results_M_5000.RData"))
M_5000 <- metric_comb %>%
  filter(n > 1, r < 100) %>%
  select(Metric, q95, n_cat) %>%
  mutate(M = 5000)

M_all <- rbind(M_500, M_1000, M_5000) %>% mutate(M = factor(M))

M_all %>%
  mutate(q95 = ifelse(q95 < 1e-4, 1e-4, q95)) %>%
  filter(Metric == "w") %>%
  group_by(n_cat, M) %>%
  summarize(med = median(q95, na.rm = T),
            quant_05 = quantile(q95, probs = 0.25, na.rm = T),
            quant_95 = quantile(q95, probs = 0.75, na.rm = T)) %>%
  ggplot() +
  geom_ribbon(aes(x = as.numeric(n_cat), ymin = quant_05, ymax = quant_95, fill = M),
              alpha = 0.15) +
  geom_line(aes(x = as.numeric(n_cat), y = med, color = M), lwd = 2) +
  theme_classic() +
  labs(x = "Sample Size (n)", 
       y = "95th Percentile",
       title = "") +
  scale_color_manual(values = brewer.pal(n = 9, name = "Set1"), name = "M") +
  scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"), name = "M") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:length(levels(metric_comb$n_cat)),
                     labels = c("(1, 10]", "(10, 100]", "(100, 1k]", 
                                "(1k, 10k]", "> 10k"),
                     expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")
ggsave(filename = "q95_M_sensitivity.png", path = output_dir,
       width = 8, height = 6, unit = "in")


################################################
### Sensitivity analysis: vary window length ###
################################################

window_7 <- all_locations %>%
  filter(Location %in% c("Brazil", "Denmark", "United States"),
         n > 1) %>%
  select(Date, Location, delay_days, w) %>%
  mutate(Window = 7)

window_3 <- fread(paste0(results_dir, "all_metric_results_3_window.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2020-11-01")) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant",
         n > 1) %>%
  arrange(Date) %>%
  select(Date, Location, delay_days, w) %>%
  mutate(Window = 3)

window_14 <- fread(paste0(results_dir, "all_metric_results_14_window.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2020-11-01")) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant",
         n > 1) %>%
  select(Date, Location, delay_days, w) %>%
  mutate(Window = 14)

window_comb <- rbind(window_7, window_3, window_14)

window_comb %>%
  filter(delay_days <= 14) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(delay_days, Location, Window) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(line_group = cumsum(is.na(w))) %>%
  ggplot(aes(x = Date, y = w, color = as.factor(Window))) +
  geom_line(aes(group = interaction(delay_days, Window, line_group)),
            alpha = 0.5, linewidth = 1) +
  geom_point(alpha = 0.7, size = 2, show.legend = T) +
  scale_color_manual(values = brewer.pal(4, name = "Dark2"),
                     name = "Window \nLength",
                     na.translate = F) +
  theme_bw() +
  labs(x = "", 
       y = "Cohen's w",
       title = "") +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  ggh4x::facet_grid2(
    Location ~ delay_days,
    scales = "free_y",
    independent = "y")

ggsave(filename = "Window_sensitivity.png", path = output_dir,
       width = 15, height = 10, unit = "in")

####################################
### Collapsed Omicron categories ###
####################################

full_omicron <- fread(paste0(results_dir, "all_metric_results.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2021-12-01")) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant",
         Location %in% c("Brazil", "Denmark", "United States")) %>%
  select(Date, Location, delay_days, w) %>%
  arrange(Date) %>%
  mutate(omicron_categories = "Full")

collapse_omicron <- fread(paste0(results_dir, "all_metric_results_collapse_omicron.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2021-12-01")) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant") %>%
  select(Date, Location, delay_days, w) %>%
  arrange(Date) %>%
  mutate(omicron_categories = "Collapsed")

omicron_comb <- rbind(full_omicron, collapse_omicron)

### time series plots ###
omicron_comb %>%
  filter(delay_days <= 14) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(Location, delay_days, omicron_categories) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(line_group = cumsum(is.na(w))) %>%
  ungroup() %>%
  ggplot(aes(x = Date, y = w, color = omicron_categories)) +
  geom_line(aes(group = interaction(Location, delay_days, omicron_categories, line_group)),
            alpha = 0.5, linewidth = 0.5) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_classic() +
  labs(x = "", y = "Cohen's w") +
  scale_color_manual(
    values = brewer.pal(4, name = "Dark2"),
    name = "Omicron \nCategories",
    na.translate = FALSE) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  ggh4x::facet_grid2(
    Location ~ delay_days,
    scales = "free_y",
    independent = "y"
  )

ggsave(filename = "category_sensitivity.png", path = output_dir,
       width = 12, height = 6, unit = "in")

### correlation ###
omicron_wide <- omicron_comb %>%
  pivot_wider(names_from = "omicron_categories",
              values_from = "w")

find_cor <- function(loc, delay, method = "pearson"){
  sub_df <- omicron_wide %>%
    filter(Location == loc,
           delay_days == delay)
  
  cor_val <- cor(sub_df$Full, sub_df$Collapsed,
                 use = 'pairwise.complete.obs',
                 method = method)
  
  return(round(cor_val, 3))
}

find_cor("Brazil", 14)
find_cor("Denmark", 7)
find_cor("Denmark", 14)
find_cor("United States", 7)
find_cor("United States", 14)

find_cor("Brazil", 7, "spearman")
find_cor("Brazil", 14, "spearman")
find_cor("Denmark", 7, "spearman")
find_cor("Denmark", 14, "spearman")
find_cor("United States", 7, "spearman")
find_cor("United States", 14, "spearman")

find_cor("Brazil", 21)
find_cor("Brazil", 30)
find_cor("Denmark", 21)
find_cor("Denmark", 30)
find_cor("United States", 21)
find_cor("United States", 30)

find_cor("Brazil", 21, "spearman")
find_cor("Brazil", 30, "spearman")
find_cor("Denmark", 21, "spearman")
find_cor("Denmark", 30, "spearman")
find_cor("United States", 21, "spearman")
find_cor("United States", 30, "spearman")

#################################
### Further Fei investigation ###
#################################

load(paste0(results_dir, "all_metric_results.RData"))

Fei_df <- for_sim %>%
  mutate(min_pval = sapply(for_sim$p_val, min)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2020-11-01"),
         n > 1,
         r < 100)

### Fei vs min p_val ###
all_min_p <- Fei_df %>%
  mutate(n_int = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("n: 2-50", "n: 50-500", "n: 500-5,000", "n > 5,000"),
                     include.lowest = T),
         K_int = cut(n,
                     breaks = c(2, 3, 4, 5, 6, Inf),
                     labels = c("2", "3", "4", "5", "> 6"),
                     include.lowest = T),
         pval_int = cut(min_pval,
                        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5),
                        labels = c("(0, 0.01]", "(0.01, 0.05]", "(0.05, 0.10]",
                                   "(0.1, 0.15]", "(0.15, 0.20]", "(0.20, 0.25]",
                                   "(0.25, 0.5]"))) %>%
  ggplot() +
  geom_jitter(aes(x = pval_int, y = Fei, color = K_int),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = pval_int, y = Fei),
               alpha = 0) +
  scale_color_viridis_d(name = "K",
                        option = "plasma") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic() +
  labs(x = "Min pval", 
       y = "Fei") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ n_int, nrow = 1)

min_p_0.01 <- Fei_df %>%
  filter(min_pval <= 0.01) %>%
  mutate(n_int = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("n: 2-50", "n: 50-500", "n: 500-5,000", "n > 5,000"),
                     include.lowest = T),
         K_int = cut(n,
                     breaks = c(2, 3, 4, 5, 6, Inf),
                     labels = c("2", "3", "4", "5", "> 6"),
                     include.lowest = T),
         pval_int = cut(min_pval,
                        breaks = c(0, 1e-4, 5e-4, 1e-3, 5e-3, 0.01),
                        labels = c("(0, 1e-4]", "(1e-4, 5e-4]", "(5e-4, 1e-3]",
                                   "(1e-3, 5e-3]", "(5e-3, 0.01]"),
                        include.lowest = T)) %>%
  ggplot() +
  geom_jitter(aes(x = pval_int, y = Fei, color = K_int),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = pval_int, y = Fei),
               alpha = 0) +
  scale_color_viridis_d(name = "K",
                        option = "plasma") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic() +
  labs(x = "Min pval", 
       y = "Fei") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ n_int, nrow = 1)

min_p_0.01_trunc <- Fei_df %>%
  filter(min_pval <= 0.01) %>%
  mutate(n_int = cut(n,
                     breaks = c(0, 50, 500, 5000, Inf),
                     labels = c("n: 2-50", "n: 50-500", "n: 500-5,000", "n > 5,000"),
                     include.lowest = T),
         K_int = cut(n,
                     breaks = c(2, 3, 4, 5, 6, Inf),
                     labels = c("2", "3", "4", "5", "> 6"),
                     include.lowest = T),
         pval_int = cut(min_pval,
                        breaks = c(0, 1e-4, 5e-4, 1e-3, 5e-3, 0.01),
                        labels = c("(0, 1e-4]", "(1e-4, 5e-4]", "(5e-4, 1e-3]",
                                   "(1e-3, 5e-3]", "(5e-3, 0.01]"),
                        include.lowest = T)) %>%
  ggplot() +
  geom_jitter(aes(x = pval_int, y = Fei, color = K_int),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = pval_int, y = Fei),
               alpha = 0) +
  scale_color_viridis_d(name = "K",
                        option = "plasma") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic() +
  labs(x = "Min pval", 
       y = "Fei") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ n_int, nrow = 1) +
  coord_cartesian(ylim = c(0, 0.1))

plot_grid(all_min_p, min_p_0.01, min_p_0.01_trunc, ncol = 1, align = "v", axis = "l")
ggsave(filename = "Fei_vs_min_pval.png", path = output_dir,
       width = 10, height = 15, unit = "in")

### Fei vs Cohen's w ###
cor(Fei_df$w, Fei_df$Fei, use = "pairwise.complete.obs")
cor(Fei_df$w, Fei_df$Fei, use = "pairwise.complete.obs",
    method = "spearman")

Fei_df %>%
  mutate(pval_int = cut(min_pval,
                        breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5),
                        labels = c("(0, 0.01]", "(0.01, 0.05]", "(0.05, 0.10]",
                                   "(0.1, 0.25]", "(0.25, 0.5]"),
                        include.lowest = T)) %>%
  ggplot() +
  geom_point(aes(x = Fei, y = w, color = min_pval),
             alpha = 0.1, size = 1.5, show.legend = T) +
  theme_classic() +
  labs(x = "Fei",
       y = "Cohen's w",
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_viridis(option = "viridis",
                      name = "Min pval",
                      limits = c(0, 0.5),
                      direction = -1) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
ggh4x::facet_grid2(
  ~ pval_int,
  scales = "free_y",
  independent = "y"
)

ggsave(filename = "w_vs_Fei_scatter.png", path = output_dir,
       width = 15, height = 4, unit = "in")
