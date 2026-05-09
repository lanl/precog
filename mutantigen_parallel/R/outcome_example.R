library(tidyverse)
library(pracma) # for finding peaks
### Data examples

data_dir <- "../mutantigen_experiment/outfiles/"
trial_num <- "1"

#### Calculate Epi Outcomes

ts <- read.delim(paste0(data_dir, "out_", trial_num, ".timeseries"))

# 1. Attack Rate (this over a rolling window would be better)
ts$ar      <- ts$totalI/ts$totalS*100
average_ar <- mean(ts$ar)
sd_ar      <- sd(ts$ar)

# 2. Number of peaks
totalI_smooth <- loess(ts$totalI ~ ts$date, span = .05)
ts$totalI_smooth <- predict(totalI_smooth)

peaks <- findpeaks(ts$totalI_smooth, minpeakdistance = 15, nups = 5, minpeakheight = 300, npeaks = 20)
peak_heights <- peaks[,1]

# code for checking heights
#peak_times <- data.frame(date = ts$date[peaks[, 2]])
#ts %>% 
#  ggplot(aes(x = date, y = totalI_smooth)) + geom_line() + 
#  geom_line(aes(x = date, y = totalI), color = "blue") +
#  geom_vline(data = peak_times, aes(xintercept = date))

npeaks <- nrow(peaks)
sd_heights <- sd(peak_heights)    

if(is.null(npeaks)) {
  npeaks <- 0
  sd_heights <- 0
}

##### Calculate Evo Outcomes
antigen_ts      <-
  read_csv(
    paste0(data_dir, "out_", trial_num, ".antigenicSamples"),
    show_col_types = FALSE,
    name_repair = "universal_quiet"
  )

# 1. Number of VOC over 20 years
num_voc      <- ncol(antigen_ts) - 1
colnames(antigen_ts)[1] <- "time"
colnames(antigen_ts)[2:ncol(antigen_ts)] <- paste0("voc_", seq(1:num_voc))


# 2. Is there turnover? 
adequate_samples <- which(rowSums(antigen_ts[, -1]) > 0) # Sometimes last time point doesn't have any sequences sampled if I is low
turnover <- as.numeric(ifelse(which.max(antigen_ts[tail(adequate_samples, 1), 2:ncol(antigen_ts)]) == 1, 0, 1))

# 3. Median VOC circulating 
antigen_ts <- antigen_ts %>% 
  pivot_longer(names_to = "voc", values_to = "numbers_sampled", -time, names_prefix = "voc_") 
antigen_ts$voc <- factor(antigen_ts$voc, levels = seq(1:num_voc))

# First - Turn into proportion and apply threhsold of 0.1
voc_prop <- antigen_ts %>% 
  group_by(time) %>% 
  mutate(rel_freq = numbers_sampled/sum(numbers_sampled)) %>% 
  select(time, voc, rel_freq) 

median_circulating <- voc_prop %>% 
  group_by(time) %>% 
  filter(rel_freq > 0.1) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  summarize(median_circulating = median(n),
            mean_circulating = mean(n)) %>% 
  select(median_circulating)

# Put together 
outcomes_profile <- data.frame(ar= average_ar, num_voc = num_voc, sd_epi_heights = sd_heights,
                          num_epi_peaks = npeaks, medidan_circulating = median_circulating$median_circulating, 
                          turnover = turnover)

