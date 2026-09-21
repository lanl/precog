#rm(list=ls())

#########################
### Main Body Figures ###
#########################

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
output_dir = paste0(dir, "plots/main_body/")

##############################
### Read in and clean data ###
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

### dominant variant/variant turnover ###
dominant_variant_w <- read_csv(
  paste0(results_dir, "dominant_variant_w_merged.csv"),
  show_col_types = FALSE
) %>%
  mutate(Date = as.Date(Date))

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



#################
### w vs n, r ###
#################

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
  filter(Metric == "w",
         r < 100,
         n > 1,
         #n <= 5000,
         True_Value >= 1,
  )
num <- denom %>% filter(pvalue < 0.05)

nrow(num)/nrow(denom) 



#################################
### Simulations & Time Series ###
#################################

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

make_emergence_plot(loc = "United States", delay = 14)
make_emergence_plot(loc = "Denmark", delay = 14)

########################
### Variant turnover ###
########################

plot_data <- dominant_variant_w %>%
  filter(!is.na(delta_p_val),
         !is.na(w),
         n > 1,
         r < 100) %>%
  mutate(delay_days = factor(delay_days, levels = c(7, 14, 21, 30)),
    n_int = cut(n,
                breaks = c(0, 50, 500, Inf),
                labels = c("(1, 50]", "(50, 500]", "> 500"),
                include.lowest = TRUE),
    r_int = cut(r,
                breaks = c(0, 10, 25, 50, 100),
                labels = c("(0, 10]", "(10, 25]", "(25, 50]", "(50, 100)"),
                include.lowest = TRUE),
    delta_int = cut(delta_p_val,
                    breaks = c(-1, -0.25, -0.1, -0.05, -0.025, -0.01,
                               0, 0.01, 0.025, 0.05, 0.1, 0.25, 1),
                    labels = c("[-1, -0.25]", "(-0.25, -0.1]", "(-0.1, -0.05]",
                               "(-0.05, -0.025]", "(-0.025, -0.01]", "(-0.01, 0]",
                               "(0, 0.01]", "(0.01, 0.025]", "(0.025, 0.05]",
                               "(0.05, 0.1]", "(0.1, 0.25]", "(0.25, 1]"),
      include.lowest = TRUE),
    abs_delta = abs(delta_p_val),
    abs_delta_int = cut(abs_delta,
                        breaks = c(0, 0.01, 0.025, 0.05, 0.1, 1),
                        labels = c("(0, 0.01]", "(0.01, 0.025]", 
                                   "(0.025, 0.05]", "(0.05, 0.1]", 
                                   "(0.1, 1]"),
                        include.lowest = TRUE))

### daily rate of change of dominant variant ###
abs_rate_w_plot <- ggplot(plot_data) +
  geom_jitter(
    aes(x = factor(abs_delta_int), y = w, color = n_int),
    alpha = 0.1,
    show.legend = TRUE) +
  geom_boxplot(aes(x = factor(abs_delta_int), y = w), alpha = 0) +
  labs(
    x = "Absolute daily rate of change in dominant variant proportion",
    y = "Cohen's w") +
  scale_color_viridis_d(
    option = "plasma",
    direction = -1,
    name = "n" ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  ) +
  facet_grid(~ delay_days, labeller = label_both) +
  coord_cartesian(ylim = c(0, 2))

ggsave(
  paste0(output_dir, "dominant_variant.png"),
  abs_rate_w_plot,
  width = 9,
  height = 4,
  dpi = 300
)

# Numerical summaries use the same row filters as the plotting script.
correlation_summary <- dominant_variant_w %>%
  filter(
    !is.na(delta_p_val),
    !is.na(w),
    n > 1,
    r < 100
  ) %>%
  group_by(delay_days) %>%
  summarise(
    n = n(),
    pearson_r = cor(abs(delta_p_val), w, method = "pearson"),
    spearman_rho = cor(abs(delta_p_val), w, method = "spearman"),
    .groups = "drop"
  )

### variant turnover ###
prepare_turnover_data <- function(data, turnover_column) {
  data %>%
    filter(delay_days == 7) %>%
    mutate(
      turnover = factor(
        .data[[turnover_column]],
        levels = c(FALSE, TRUE),
        labels = c("Non-turnover", "Turnover")),
      n_int = cut(n,
                  breaks = c(0, 50, 500, Inf),
                  labels = c("n: 2-50", "n: 51-500", "n > 500"),
                  include.lowest = TRUE)) %>%
    group_by(n_int, r_int, turnover) %>%
    mutate(
      q1 = quantile(w, 0.25, na.rm = TRUE),
      q3 = quantile(w, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      plot_point = between(w, q1 - 1.5 * iqr, q3 + 1.5 * iqr)
    ) %>%
    ungroup()
}

turnover_cols <- c(
  "Non-turnover" = "#56B4E9",
  "Turnover"     = "#E69F00"
)

plot_turnover <- function(turnover_data, fig_name) {
  p <- ggplot(turnover_data, aes(x = turnover, y = w)) +
    geom_jitter(
      data = ~ filter(.x, plot_point),
      aes(color = turnover),
      alpha = 0.15
    ) +
    geom_boxplot(aes(fill = turnover), outliers = FALSE, alpha = 0) +
    scale_color_manual(values = turnover_cols) +
    scale_fill_manual(values = turnover_cols) +
    ggh4x::facet_grid2(
      ~ n_int,
      scales = "free_y",
      independent = "y"
    ) +
    labs(x = NULL, y = "Cohen's w") +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(
    paste0(output_dir, fig_name),
    p,
    width = 9,
    height = 3,
    dpi = 300
  )
}

turnover_2wk_before_after <- prepare_turnover_data(
  plot_data,
  "turnover_2wk_before_after"
)

plot_turnover(turnover_2wk_before_after, "variant_turnover.png")
