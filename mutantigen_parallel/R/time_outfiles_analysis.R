# Look through time output files of Mutantigen Parallel to analyze the runtimes
# of the sims that do complete.
## Author: AC Murph
## Date: February 2025
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(tidyverse)
library(this.path)
library(stringr)
library(readr)
library(latex2exp)
library(cowplot)
library(png)
library(ggtext)
library(sysfonts)
library(showtext)
setwd(this.path::here())

# Iterate through all time_logfiles expected
all_outfiles = list.files('../time_logfiles/')
input_csv = as_tibble(read.csv("../lhs_samples.csv"))
all_outfiles_split = strsplit(all_outfiles, split = '\\.')
scale = 1.4
has_all_files = c()
has_timeseries_not_trees = c()
num_of_files = c()
full_file_df = NULL
all_numbers = c()
all_large_pops = c()
all_attack_rates = c()
all_file_exists = c()


number_of_reps = 20

good_input_files = as_tibble(read.csv('../antipeak.csv'))

for(curr_file in all_outfiles){
  # Define the file path
  file_path <- paste0('../time_logfiles/', curr_file)
  
  i <- str_extract(curr_file, "(?<=_)\\d+(?=\\.)") |> as.integer()
  
  num_of_yaml_row = (i-1) %% nrow(good_input_files) + 1
  num_of_yaml = good_input_files$runnum[num_of_yaml_row]
  
  name_of_yaml = paste0("input_files/parameters_load_updated_", num_of_yaml, ".yml")
  
  # Read the third line
  lines <- readLines(file_path)
  
  # also read in timeseries file.
  path <- paste0("../outfiles/out_",i, ".timeseries")  # <- change to your filename/path
  second_path = paste0("../outfiles/out_",i,".out.simmapAntigenic")
  
  # Ensure there are at least three lines
  if (length(lines) >= 3) {
    third_line <- lines[3]  # R uses 1-based indexing
    
    # Use a regular expression to extract the number after ":"
    match <- regmatches(third_line, regexpr(":(\\s*-?\\d+\\.?\\d*)", third_line))
    
    if (length(match) > 0) {
      # Extract only the numeric part and convert it to numeric
      number <- as.numeric(gsub(":", "", match))
      print(paste("Extracted number:", number))
      all_numbers = c(all_numbers, number)
      xx = input_csv$initialNs[num_of_yaml]>7.5e6
      if(is.na(xx)) browser()
      all_large_pops = c(all_large_pops, xx)
      
      if(file.exists(path)){
        # Read as TSV; default everything to numeric (double)
        df <- read_tsv(
          file = path,
          col_types = cols(.default = col_double()),
          progress = FALSE,
          trim_ws = TRUE
        )
        max_attack_rate = max(df$totalI/(df$totalS+df$totalI))
        all_attack_rates = c(all_attack_rates, max_attack_rate)
      } else {
        all_attack_rates = c(all_attack_rates, NA)
      }
      
      
      all_file_exists = c(all_file_exists, (file.exists(path))&(file.exists(second_path)))
    } else {
      print("No number found after ':' in the third line.")
    }
  } else {
    print("File has less than three lines.")
  }
  
}

df <- data.frame(x = all_numbers, large_population = all_large_pops) %>%
  mutate(
    large_population = factor(
      large_population,
      levels = c(TRUE, FALSE),
      labels = c("(> 7.5 million)", "(<= 7.5 million)")
    )
  )

# Helper to render LaTeX in the legend/axis
labs_tex <- function(vals) {
  map <- c(
    "(> 7.5 million)"  = "$> 7.5 million$",
    "(<= 7.5 million)" = "$\\leq 7.5 million$"
  )
  latex2exp::TeX(map[vals])
}

# 2) A more visually appealing ggplot (adjust x/y & geoms to your data)
p <- ggplot(df, aes(x = x, fill = large_population)) +
  geom_histogram(linewidth = 2, alpha = 0.85, bins = 50) +
  # tasteful palette + LaTeX legend labels
  scale_fill_manual(
    values = c("(> 7.5 million)" = "#3366CC",
               "(<= 7.5 million)" = "#DC3912"),
    labels = labs_tex,
    name = NULL
  ) +
  labs(
    title = "Runtimes by InitialNs",
    subtitle = latex2exp::TeX("Colored by $> 7.5$ million vs. $\\leq 7.5$ million"),
    x = "Time (in minutes)",
    y = "Counts"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 6))
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))+
  annotate("text", x = max(all_numbers) * 0.7, y = 1500, 
           label = paste("Number of sims completed \nwithout error:", length(all_numbers)), size = 5, color = "steelblue")+
  # BIGGER BASE SIZE
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    
    # Legend emphases
    legend.position = "top",
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.key.size = unit(10, "mm"),
    
    # Titles
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18, margin = margin(b = 6))
  ) +
  # make legend points bigger too
  guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1))) + 
  geom_vline(xintercept = 8*60, linetype = 'dashed')

png("time_analysis_InitialNs.png", width = 700 * scale, height = 480 * scale)
p
dev.off()


df <- data.frame(x = all_numbers, 
                 max_attack_rate = all_attack_rates, 
                 both_outfiles_exist = all_file_exists)
# OPTIONAL: path to a hype background image (speed lines, burst, etc.)
# e.g., "speedlines.png" (1920x1080 with transparency looks great)
bg_path <- "speedlines.png"  # <- replace with your file; leave as-is to skip

title_html <- paste0(
  "Runtimes by Max ",
  "<span style='",
  "font-family: Impact, Haettenschweiler, sans-serif;",
  "font-weight: 900;",
  "letter-spacing: 1px;",
  "color: #ff1744;",
    "text-shadow:",
  " 0 0 1px #ffffff,",   # thin white halo
  " 0 0 6px #ff5252,",   # red glow
  " 0 0 14px #ff9100,",  # orange glow
  " -3px -3px 0 #000, 3px -3px 0 #000, -3px 3px 0 #000, 3px 3px 0 #000,", # outline
  " -5px 0 0 #000, 5px 0 0 #000, 0 -5px 0 #000, 0 5px 0 #000;",
  "transform: skewX(-6deg);",
  "'>ATTACK</span>    Rate"
)

p_base <-
  ggplot(df, aes(x = x, y = max_attack_rate)) +
  geom_point(aes(color = both_outfiles_exist)) +
  labs(
    title = title_html,
    subtitle = latex2exp::TeX("Determined as max( I / (S + I) )"),
    x = "Time (in minutes)",
    y = "Max Attack Rate"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = ggtext::element_markdown(
      size = 22,                   # base size (inline CSS can push bigger)
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(size = 18, margin = margin(b = 6)),
    text = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.key.size = unit(10, "mm")
  ) +
  guides(color = guide_legend(title = "Both outfiles exist:",
                              override.aes = list(size = 4, alpha = 1))) +
  annotate(
    "text",
    x = max(all_numbers) * 0.7,
    y = max(all_attack_rates) - min(all_attack_rates),
    label = paste("Number of sims completed \nwithout error:", length(all_numbers)),
    size = 5,
    color = "steelblue"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  geom_vline(xintercept = 8 * 60, linetype = "dashed")+ 
  theme(
    plot.background  = element_rect(fill = NA, color = NA),
    panel.background = element_rect(fill = NA, color = NA)
  )

# If you provided a background image, composit it under the plot.
# Transparent PNGs with radial burst or speedlines work best.
final_plot = p_base
final_plot

png("time_analysis_MaxAttackRates.png", width = 700 * scale, height = 480 * scale)
final_plot
dev.off()

