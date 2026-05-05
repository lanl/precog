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

dir = "./reporting-delay-pub/"

##############################
### Read in data and clean ###
##############################

### GISAID data for data section figure ###
merged_df <- read.csv(paste0(dir, "reporting_delay_data.csv")) %>%
  select(-X) %>%
  rename(Date = "collection_date",
         Location = "Admin0") %>%
  mutate(Year = year(Date)) %>%
  filter(Location %in% c("Global", "United States", "Denmark", "Brazil")) %>%
  mutate(Location = factor(Location, levels = c("Brazil", "Denmark", "United States", "Global")))

### metric results for all locations ###
all_locations <- read.csv(file = paste0(dir, "all_cramer_results.csv")) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date < as.Date("2023-01-01"),
         Date >= as.Date("2020-11-01")) %>%
  mutate(error_code = ifelse(is.na(error_code), "0", error_code)) %>%
  filter(error_code != "All samples are of the same variant") %>%
  select(-c(X, error_code)) %>%
  arrange(Date) %>%
  mutate(L1 = L1/2,
         L2 = L2/sqrt(2),
         omega_cat = cut(omega, 
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
load(file = paste0(dir, "all_sim_results.RData")) # this is metric_comb


### results for emerging variant time series ###
emerge_join <- read.csv(file = paste0(dir, "emerging_variant_ts.csv")) %>%
  mutate(Date = as.Date(collection_date)) %>%
  rename(Location = Admin0,
         delay_days = delay_bin) %>%
  select(-c(X, collection_date)) %>%
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
ggsave(filename = "compare_delays.png", path = dir,
       width = 12, height = 4, unit = "in", limitsize = F)

### percent reported within 1 and 2 weeks ###

# disaggregated by year
merged_df %>%
  filter(Year < 2023) %>%
  group_by(Admin0, Year) %>%
  summarize(n_total = sum(counts, na.rm = T),
            n_1wk = sum(counts[delay_days <= 7], na.rm = T),
            pct_1wk = 100 * (n_1wk/n_total),
            n_2wk = sum(counts[delay_days <= 14], na.rm = T),
            pct_2wk = 100 * (n_2wk/n_total),
            .groups = "drop") %>%
  select(-n_1wk, -n_2wk)

merged_df %>%
  filter(Year < 2023) %>%
  group_by(Admin0) %>%
  summarize(n_total = sum(counts, na.rm = T),
            n_1wk = sum(counts[delay_days <= 7], na.rm = T),
            pct_1wk = 100 * (n_1wk/n_total),
            n_2wk = sum(counts[delay_days <= 14], na.rm = T),
            pct_2wk = 100 * (n_2wk/n_total),
            .groups = "drop") %>%
  select(-n_1wk, -n_2wk)

##########################
### n-omega categories ###
##########################

### all locations, all delay periods ###

#categories align with CV_vs_reporting_rate.png
n_omega_cats <- metric_comb %>% 
  filter(n > 1, omega < 100) %>%
  mutate(n_cat = cut(n,
                     breaks = c(0, 100, 1000, 10000, Inf),
                     labels = c("[2, 100]", "(100, 1k]", 
                                "(1k, 10k]", "> 10k"), 
                     include.lowest = T)) %>%
  count(n_cat, omega_cat) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = omega_cat, y = n_cat, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1),
                color = prop > 0.15), size = 5) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "Sample Size (n)") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Proportion of Observations",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")


# statistial significance
n_omega_sig <- metric_comb %>% 
  filter(n > 1, omega < 100,
         Metric == "V") %>%
  mutate(n_cat = cut(n,
                     breaks = c(0, 100, 1000, 10000, Inf),
                     labels = c("[2, 100]", "(100, 1k]", 
                                "(1k, 10k]", "> 10k"), 
                     include.lowest = T)) %>%
  count(n_cat, omega_cat, pvalue_cat) %>%
  pivot_wider(names_from = pvalue_cat, values_from = n) %>%
  mutate(sig_prop = `< 0.05` / (`< 0.05` + `> 0.05`),
         sig_prop = case_when(is.na(`< 0.05`) ~ 0,
                              is.na(`> 0.05`) ~ 1,
                              .default = sig_prop)) %>% 
  ggplot(aes(x = omega_cat, y = n_cat, fill = sig_prop)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(sig_prop, accuracy = 0.1),
                color = sig_prop > 0.2), size = 5) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Significance Proportion",
                       na.value = "white", option = "magma",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")


plot_grid(n_omega_cats, n_omega_sig, ncol = 2, align = "h", axis = "l")
ggsave(filename = "n_omega_cats_all_countries.png", path = dir,
       width = 16, height = 6, unit = "in")

### separate delay periods ###

all_countries_cat <- metric_comb %>% 
  filter(n > 1, #omega < 100, 
         Location != "Global",
         Metric == "V") %>%
  mutate(omega_cat = cut(omega,
                    breaks = c(0, 20, 40, 60, 80, 100),
                    labels = c("(0, 20]", "(20, 40]", "(40, 60]",
                               "(60, 80]", "(80, 100]"),
                    include.lowest = T),
         n_cat = cut(n,
                     breaks = c(0, 10, 100, 1000, 10000, Inf),
                     labels = c("[2, 10]", "(10, 100]", "(100, 1k]", "(1k, 10k]",
                                "> 10k"),
                     include.lowest = T)) %>%
  group_by(delay_days) %>%
  count(n_cat, omega_cat) %>%
  mutate(prop = n/ sum(n),
         Location = "All Countries") %>%
  select(Location, delay_days, n_cat, omega_cat, n, prop)


by_location_cat <- metric_comb %>% 
  filter(n > 0, #omega < 100, 
         Metric == "V",
         Location %in% c("Brazil", "Denmark", "United States")) %>%
  mutate(omega_cat = cut(omega,
                         breaks = c(0, 20, 40, 60, 80, 100),
                         labels = c("(0, 20]", "(20, 40]", "(40, 60]",
                                    "(60, 80]", "(80, 100]"),
                         include.lowest = T),
         n_cat = cut(n,
                     breaks = c(0, 10, 100, 1000, 10000, Inf),
                     labels = c("[2, 10]", "(10, 100]", "(100, 1k]", 
                               "(1k, 10k]","> 10k"),
                     include.lowest = T)) %>%
  group_by(Location, delay_days) %>%
  count(n_cat, omega_cat) %>%
  mutate(prop = n/ sum(n))
  
# combine for plot
rbind(all_countries_cat, by_location_cat) %>%
  ggplot(aes(x = omega_cat, y = n_cat, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1),
                color = prop > 0.3), size = 3) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "Sample Size (n)") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Proportion of Observations",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_grid(Location ~ delay_days)
ggsave(filename = "n_omega_cats.png", path = dir,
       width = 10, height = 10, unit = "in")


# statistical significance
all_countries_sig <- metric_comb %>% 
  filter(n > 1, omega < 100, 
         Location != "Global",
         Metric == "V") %>%
  mutate(omega_cat = cut(omega,
                         breaks = c(0, 20, 40, 60, 80, 100),
                         labels = c("(0, 20]", "(20, 40]", "(40, 60]",
                                    "(60, 80]", "(80, 100]"),
                         include.lowest = T),
         n_cat = cut(n,
                     breaks = c(0, 10, 100, 1000, 10000, Inf),
                     labels = c("[2, 10]", "(10, 100]", "(100, 1k]", 
                                "(1k, 10k]","> 10k"),
                     include.lowest = T)) %>%
  group_by(delay_days) %>%
  count(n_cat, omega_cat, pvalue_cat) %>%
  pivot_wider(names_from = pvalue_cat, values_from = n) %>%
  mutate(sig_prop = `< 0.05` / (`< 0.05` + `> 0.05`),
         sig_prop = case_when(is.na(`< 0.05`) ~ 0,
                              is.na(`> 0.05`) ~ 1,
                              .default = sig_prop)) %>% 
  mutate(Location = "All Countries") %>%
  select(Location, delay_days, n_cat, omega_cat, `< 0.05`, `> 0.05`, sig_prop)


by_location_sig <- metric_comb %>% 
  filter(n > 1, omega < 100, Metric == "V",
         Location %in% c("Brazil", "Denmark", "United States")) %>%
  mutate(omega_cat = cut(omega,
                         breaks = c(0, 20, 40, 60, 80, 100),
                         labels = c("(0, 20]", "(20, 40]", "(40, 60]",
                                    "(60, 80]", "(80, 100]"),
                         include.lowest = T),
         n_cat = cut(n,
                     breaks = c(0, 10, 100, 1000, 10000, Inf),
                     labels = c("[2, 10]", "(10, 100]", "(100, 1k]", 
                                "(1k, 10k]","> 10k"),
                     include.lowest = T)) %>%
  group_by(Location, delay_days) %>%
  count(n_cat, omega_cat, pvalue_cat) %>%
  pivot_wider(names_from = pvalue_cat, values_from = n) %>%
  mutate(sig_prop = `< 0.05` / (`< 0.05` + `> 0.05`),
         sig_prop = case_when(is.na(`< 0.05`) ~ 0,
                              is.na(`> 0.05`) ~ 1,
                              .default = sig_prop)) 


rbind(all_countries_sig, by_location_sig) %>%
  ggplot(aes(x = omega_cat, y = n_cat, fill = sig_prop)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(sig_prop, accuracy = 0.1),
                color = sig_prop > 0.2), size = 3) +
  theme_bw() +
  labs(x = "Reporting Rate (%)", y = "") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis_c(direction = -1, name = "Significance Proportion",
                       na.value = "white", option = "magma",
                       guide = guide_colorbar(barwidth = 15, barheight = 1)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_grid(Location ~ delay_days)
ggsave(filename = "n_omega_sig.png", path = dir,
       width = 10, height = 10, unit = "in")

###########################  
### Cramer's V vs Norms ###
###########################

## when considering observations with at least 2 near-real-time samples
filter_n <- all_locations %>% filter(n > 1)

cor(filter_n$V, filter_n$L1,
    use = 'pairwise.complete.obs') %>% round(3)
cor(filter_n$V, filter_n$L2,
    use = 'pairwise.complete.obs') %>% round(3)
cor(filter_n$V, filter_n$Linf,
    use = 'pairwise.complete.obs') %>% round(3)

### boxplot ###
p <- all_locations %>%
  filter(n > 1,
         !(is.na(V))) %>%
  select(V, L1, L2, Linf, omega) %>%
  pivot_longer(cols = c(L1, L2, Linf), names_to = "Norm", values_to = "Norm_Value") %>%
  mutate(norm_int = cut(Norm_Value, 
                        breaks = seq(0, 1, 0.25), 
                        labels = c("(0, 0.25]", "(0.25, 0.5]", "(0.5, 0.75]", "(0.75, 1]"),
                        include.lowest = T)
  ) %>%
  ggplot(aes(x = norm_int, y = V, color = omega)) +
  geom_point(position = position_jitterdodge(), alpha = 0.2, size = 1, show.legend = T) +
  geom_boxplot(alpha = 0) +
  theme_classic() +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "Reporting Rate (%)") +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = 20, barheight = 1)) +
  labs(x = "Scaled Norm", 
       y = "Cramer\'s V") +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "bottom") +
  ylim(c(0, 1)) +
  facet_wrap(~ Norm, ncol = 1)
ggsave(filename = "Cramer_vs_Norms_box.png", p,
       path = dir, width = 7, height = 8, unit = "in")


#################################  
### Counter-intuitive results ###
#################################

cv_0 <- all_locations %>%
  filter(V == 0,
         L1 >= 0.25,
         n > 1) %>%
  select(V, n, omega, L1) %>%
  ggplot(aes(x = n, y = omega, color = L1)) +
  geom_point(alpha = 0.25, size = 3, show.legend = F) + 
  theme_classic() +
  scale_color_viridis(option = "viridis",
                      name = "Scaled L1",
                      limits = c(0, 1)) +
  #ylim(c(0, 100)) +
  labs(x = "Sample Size (n)", 
       y = "Reporting Rate (%)",
       title = "Cramer\'s V = 0 and L1 > 0.25") +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

cv_1 <- all_locations %>%
  filter(V == 1,
         L1 <= 0.25,
         n > 1) %>%
  select(V, n, omega, L1) %>%
  mutate(n_int = cut(n,
                     breaks = c(0, 10, 100, 1000),
                     labels = c("[2, 10]", "(10, 100]", "(100, 1,000]"),
                     include.lowest = T)) %>%
  ggplot(aes(x = n_int, y = omega, color = L1)) +
  geom_point(position = position_jitterdodge(), alpha = 0.1, size = 2, show.legend = T) +
  geom_boxplot(alpha = 0) +
  theme_classic() +
  scale_color_viridis(option = "viridis",
                      name = "Scaled L1",
                      limits = c(0, 1)) +
  #ylim(c(0, 100)) +
  labs(x = "Sample Size (n)", 
       y = "Reporting Rate (%)",
       title = "Cramer\'s V = 1 and L1 < 0.25") +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

grob <- arrangeGrob(cv_0, cv_1, nrow = 1)
ggsave(paste0(dir, "cv_extremes.png"),
       grob, width = 10, height = 5)

###########################
### Metrics vs n, omega ###
###########################

### Cramer's V ###
metric_comb %>%
  filter(Metric == "V",
         omega < 100,
         n > 1) %>%
  mutate(n_int = cut(n,
                     breaks = c(0, 100, 1000, 10000, Inf),
                     labels = c("n: 2-100", "n: 101-1,000", "n: 1,001-10,000", "n > 10,000"),
                              include.lowest = T)) %>%
  ggplot() +
  geom_jitter(aes(x = factor(omega_cat), y = True_Value, color = pvalue_cat),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(omega_cat), y = True_Value),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "Cramer's V") +
  scale_color_manual(values = brewer.pal(n = 4, name = "Set1"),
                     name = "p-value") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ n_int, nrow = 1)
ggsave(filename = "CV_vs_reporting_rate.png", path = dir,
       width = 15, height = 5, unit = "in")

### L1 ###
all_locations %>%
  filter(n > 1) %>%
  ggplot() +
  geom_jitter(aes(x = factor(omega_cat), y = L1, color = V),
              alpha = 0.1, show.legend = T) +
  geom_boxplot(aes(x = factor(omega_cat), y = L1),
               alpha = 0) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "Scaled L1") +
  scale_color_viridis(option = "viridis",
                      direction = -1,
                      name = "V",
                      limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ n_int, nrow = 1)
ggsave(filename = "L1_vs_reporting_rate.png", path = dir,
       width = 15, height = 5, unit = "in")

####################################
### Simulations - p-values & q95 ###
####################################

### q95 vs n ###
q95_vs_n <- metric_comb %>%
  filter(n > 1, omega < 100) %>%
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
       y = "95% Percentile",
       title = "") +
  scale_color_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:length(levels(metric_comb$n_cat)),
                     labels = c("[2, 10]", "(10, 100]", "(100, 1k]", 
                                "(1k, 10k]", "> 10k"),
                     expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom")


### q95 vs omega ###
q95_vs_omega <- metric_comb %>%
  filter(n > 1, omega < 100) %>%
  mutate(q95 = ifelse(q95 < 1e-4, 1e-4, q95),
         n_cat = cut(n,
                     breaks = c(0, 100, 1000, 10000, Inf),
                     labels = c("[2, 100]", "(100, 1k]", 
                                "(1k, 10k]", "> 10k"), 
                     include.lowest = T)) %>%
  group_by(omega_cat, Metric) %>%
  summarize(med = median(q95, na.rm = T),
            quant_05 = quantile(q95, probs = 0.25, na.rm = T),
            quant_95 = quantile(q95, probs = 0.75, na.rm = T)) %>%
  ggplot() +
  geom_ribbon(aes(x = as.numeric(omega_cat), ymin = quant_05, ymax = quant_95, fill = Metric),
              alpha = 0.15) +
  geom_line(aes(x = as.numeric(omega_cat), y = med, color = Metric), lwd = 2) +
  theme_classic() +
  labs(x = "Reporting Rate (%)", 
       y = "95% Percentile",
       title = "") +
  scale_color_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"), name = "Metric") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:length(levels(metric_comb$omega_cat)),
                     labels = levels(metric_comb$omega_cat),
                     expand = c(0, 0)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.spacing = unit(2, "cm"))


plot_grid(q95_vs_n, q95_vs_omega, ncol = 2, align = "hv", axis = "lr")
ggsave(filename = "q95.png", path = dir,
       width = 12, height = 6, unit = "in")

########################################
### Simulations - Metric Time Series ###
########################################

# true metric time series, alpha by significance
compare_V <- metric_comb %>%
  filter(Metric == "V",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = True_Value)) +
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
       y = "Cramer's V",
       title = "") +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  facet_grid(~ delay_days)


# n time series
compare_n <- metric_comb %>%
  filter(Metric == "V",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = n)) +
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


# omega time series
compare_omega <- metric_comb %>%
  filter(Metric == "V",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = omega)) +
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


# V - q95 time series
V_minus_q95 <- metric_comb %>%
  filter(Metric == "V",
         Location %in% c("Brazil", "Denmark", "United States"),
         delay_days %in% c(7, 14),
         n > 1) %>%
  mutate(delay_days = ifelse(delay_days == 7, "7 Days Post Collection", "14 Days Post Collection"),
         delay_days = factor(delay_days, levels = c("7 Days Post Collection",
                                                    "14 Days Post Collection"))) %>%
  mutate(dif = True_Value - q95) %>%
  group_by(Location, delay_days) %>%
  complete(Date = seq(min(Date), max(Date), by = 1)) %>%
  ggplot(aes(x = Date, y = dif)) +
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
       y = "V - 95% Percentile",
       title = "") +
  ylim(c(-0.1, 1)) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  facet_wrap(~ delay_days, ncol = 2)


plot_grid(compare_V, compare_n, compare_omega, V_minus_q95, ncol = 1, align = "v", axis = "l")
ggsave(filename = "Compare_Cramer_Locations.png", path = dir,
       width = 16, height = 16, unit = "in")

#######################################
### Emerging Variants - time series ###
#######################################

#loc = "United States"
#delay = 7

plot_emergence <- function(loc, delay, min_lim = -1, max_lim = 1, 
                           first_date = '2020-11-01',
                           last_date = '2022-12-31'){
  
  loc_title = ifelse(loc == "United States", "US", loc)
  first_date = as.Date(first_date)
  last_date = as.Date(last_date)
  
  pango_colors = c(brewer.pal(n = 8, name = "Dark2"), 
                   brewer.pal(n = 11, name = "RdYlBu")[c(2, 4, 8, 11)],
                   brewer.pal(n = 11, name = "PRGn")[c(2, 4, 8, 11)])
  
  sub_df <- emerge_join %>%
    filter(Location == loc,
           delay_days == delay) %>%
    arrange(Date)
  
  # determine variant ordering based on peak timing
  peak_max <- sub_df %>%
    group_by(pango) %>%
    slice_max(order_by = ret_prop, n = 1, with_ties = F) %>%
    select(pango, Date)
  
  max_props <- sub_df %>%
    group_by(pango) %>%
    summarize(max_ret_prop = max(ret_prop, na.rm = T),
              max_pro_prop = max(pro_prop, na.rm = T),
              overall_prop_dif = max_ret_prop - max_pro_prop) %>%
    mutate(Location = loc) %>%
    filter(max_ret_prop >= 0.05 | abs(overall_prop_dif) > 0.05) %>%
    merge(peak_max) %>%
    mutate(peak_time = Date) %>%
    arrange(peak_time)
  
  # only consider variants in max_props
  sub_df <- sub_df %>% filter(pango %in% max_props$pango)
  
  # variant order
  pango_order <- unique(c("Other/Unknown", unique(max_props$pango)))
  
  ### Cramer's V plot ###
  cramer <- sub_df %>%
    filter(Metric == "V",
           n > 2) %>%
    filter(Date %in% seq(first_date, last_date, by = 'days')) %>%
    mutate(pvalue_cat = factor(pvalue_cat, levels = c("> 0.05", "< 0.05"))) %>%
    ggplot() +
    geom_point(aes(x = Date, y = True_Value, color = pvalue_cat),
               alpha = 0.7, size = 3, show.legend = T) +
    theme_classic() +
    labs(x = "",
         y = "Cramer's V",
         title = paste0(loc, ": d = ", delay, " Days Post-Collection")) +
    scale_color_manual(values = brewer.pal(n = 3, name = "Paired"),
                       name = "p-value") +
    scale_x_date(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0.0)) +
    theme(text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14))
  
  ### normalized retrospective proportion ###
  ret_prop <- sub_df %>%
    mutate(pango = factor(pango, levels = pango_order)) %>%
    filter(Date %in% seq(first_date, last_date, by = 'days')) %>%
    ggplot() +
    geom_line(aes(x = Date, y = ret_prop, color = pango), 
              lwd = 1) +
    geom_point(aes(x = Date, y = ret_prop, color = pango,
                   fill = pango), alpha = 0.7, show.legend = F) +
    scale_color_manual(values = pango_colors, name = "") +
    labs(x = "", 
         y = "Validation Proportion",
         title = TeX("Validation Proportion (${p}_{val}$)")
    ) +
    theme_classic() +
    scale_x_date(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.05)) +
    theme(text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.spacing = unit(5, "cm"))
  
  ### proportion difference ###
  prop_dif <- sub_df %>%
    filter(Metric == "V") %>%
    mutate(pango = factor(pango, levels = pango_order)) %>%
    filter(Date %in% seq(first_date, last_date, by = 'days')) %>%
    ggplot() +
    geom_line(aes(x = Date, y = prop_dif, color = pango), lwd = 1) +
    geom_point(aes(x = Date, y = prop_dif, color = pango), alpha = 0.7) +
    geom_hline(yintercept = 0, color = "black", lwd = 1) +
    scale_color_manual(values = pango_colors, name = "") +
    labs(x = "Collection Date", 
         y = "Proportion Difference",
         title = TeX("Proportion Difference: ${p}_{val} - \\hat{p}_{nrt}(d)$")
    ) +
    theme_classic() +
    scale_x_date(expand = c(0, 0)) +
    theme(text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.spacing = unit(2, "cm")) +
    ylim(c(min_lim, max_lim))
  
  
  plot_grid(cramer, ret_prop, prop_dif, ncol = 1, align = "hv", axis = "lr")
  ggsave(filename = paste(loc_title, delay, "emerge.png", sep = "_"), path = dir,
         width = 12, height = 12, unit = "in")
  
}

plot_emergence("Global", 7, -0.4, 0.4)
plot_emergence("Global", 14, -0.5, 0.5)
plot_emergence("United States", 7, -0.3, 0.3)
plot_emergence("United States", 14, -0.25, 0.25)
plot_emergence("Denmark", 7, -0.5, 0.5) # remove lines for prop dif
plot_emergence("Denmark", 14, -0.2, 0.2) # remove lines for prop dif, make expand 0 for y-axis on cramer's V
plot_emergence("Brazil", 7) # remove lines for prop dif
plot_emergence("Brazil", 14) # remove lines for prop dif

