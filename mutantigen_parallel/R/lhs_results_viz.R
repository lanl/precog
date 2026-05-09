### Script to help determine parameter regimes
library(yaml)
library(tidyverse)
library(phylotate)
library(ggtree)
library(cowplot)
library(RColorBrewer)
library(pracma)

get_legend_35 <- function(plot) {
  # return all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  nonzero <- vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  # return first non-zero legend if exists, and otherwise first element (which will be a zeroGrob) 
  if (length(idx) > 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

param_ranges <- list(
  initialNs = c(10^4, 10^8),
  demeAmplitudes = c(0, 0.25),
  lambdaAntigenic = c(8.57143e-05, 0.002571429),
  meanAntigenicSize = c(1.2e-3, 1.2e-1),
  lambda = c(0.0095, 6.08),
  mutCost = c(8e-4, 8e-2),
  beta = c(0.1428, 2.25),
  nu = c(1/4, 1/14),
  epsilon_mut = c(0.5, 1.5),
  initialI_prop = c(0.0001, .001)
)



#### Step 1: Get canonical "flu" parameters:  ####


yaml_file <- "../mutantigen_parallel/parameters_load.yml"
# Load the YAML content into a list
yaml_content <- yaml.load_file(yaml_file)
  
# Convert the YAML content to a flat named vector (parameter names and values)
flat_content <- unlist(yaml_content, recursive = TRUE)
  
# Convert the vector to a data frame with parameter names as column headers
parameter_df <- as.data.frame(t(flat_content), stringsAsFactors = FALSE)
  
# Add a column for the file name
parameter_df$file <- basename(yaml_file)
flu_parameter <- parameter_df  
flu_parameter$category <- "1"

#### Step 2: Calculate flu outcomes ####
flu_profile <- data.frame()
trials <- list.dirs("../MutAntiGen/data", full.names = F)

flu_numbers <- as.numeric(str_extract(trials[grep(trials, pattern = "^trial")], "[0-9]{1,2}"))
flu_profile_df <- data.frame()
for(i in flu_numbers) {
  # Start looking at the time series and plooting what they look like
  
  dir <- paste0("../MutAntiGen/data/trial", i)
  ts_name <- paste0(dir, "/out.timeseries")
  
  ts <- read.delim(ts_name)
  ts$ar <- ts$totalI/ts$totalS*100
  average_ar <- mean(ts$ar)
  sd_ar <- sd(ts$ar)
  
  antigen_ts <- read_csv(paste0(dir, "/out.antigenicSamples"), col_names = FALSE)
  num_voc <- ncol(antigen_ts)-1
  turnover <- as.numeric(ifelse(which.max(antigen_ts[nrow(antigen_ts), 2:ncol(antigen_ts)]) == 1, 0, 1))
  
  # Number of peaks
  # Need to smooth 
  totalI_smooth <- loess(ts$totalI ~ ts$date, span = .05)
  ts$totalI_smooth <- predict(totalI_smooth)
  
  peaks <- findpeaks(ts$totalI_smooth, minpeakdistance = 15, nups = 5,  minpeakheight = 300, npeaks = 20)
  peak_heights <- peaks[,1]
  
  # code for checking heights
  peak_times <- data.frame(date = ts$date[peaks[, 2]])
  ts %>% 
    ggplot(aes(x = date, y = totalI_smooth)) + geom_line() + 
    geom_line(aes(x = date, y = totalI), color = "blue") +
    geom_vline(data = peak_times, aes(xintercept = date))
  
  npeaks <- nrow(peaks)
  sd_heights <- sd(peak_heights)    
  
  if(is.null(npeaks)) {
    npeaks <- 0
    sd_heights <- 0
  }
  flu_profile <- data.frame(ar= average_ar, num_voc = num_voc, sd_epi_heights = sd_heights, num_epi_peaks = npeaks, turnover = turnover)
  flu_profile_df <- rbind(flu_profile_df, flu_profile)
}

flu_profile_df <- cbind(flu_profile_df, flu_parameter)
flu_profile_df$disease <- "flu"
flu_outcomes <- flu_profile_df


####  Step 3: Now get all smaller set of synthetic runs ####

yaml_dir <- "../mutantigen_experiment/input_files/"

# Get a list of all YAML files in the directory
yaml_files <- list.files(yaml_dir, pattern = "\\.yml$", full.names = TRUE)

# Initialize an empty list to store data frames
parameter_list <- list()

# Loop through each YAML file
for (yaml_file in yaml_files) {
  # Load the YAML content into a list
  yaml_content <- yaml.load_file(yaml_file)
  
  # Convert the YAML content to a flat named vector (parameter names and values)
  flat_content <- unlist(yaml_content, recursive = TRUE)
  
  # Convert the vector to a data frame with parameter names as column headers
  parameter_df <- as.data.frame(t(flat_content), stringsAsFactors = FALSE)
  
  # Add a column for the file name
  parameter_df$file <- basename(yaml_file)
  
  # Append the data frame to the list
  parameter_list[[length(parameter_list) + 1]] <- parameter_df
}

# Combine all parameter data frames into a single data frame
small_parameter_df <- bind_rows(parameter_list)
small_parameter_df <- small_parameter_df %>% mutate(trial_num = as.numeric(str_extract(file,"[0-9]{1,5}")))

output_dir <- "../mutantigen_experiment/outfiles/"
figure_dir <- "../mutantigen_experiment/figures/"

# Look at the ones that have ran 
trials_success <- sort(as.numeric(str_extract(list.files(figure_dir), pattern = "[0-9]{1,4}")))
trials_error   <- seq(1:100)[-trials_success]

# What was filtered out or had an error
# What didn't run/
small_parameter_df <- small_parameter_df %>% 
  mutate(category = case_when(trial_num %in% trials_success ~ "1",
                              trial_num %in% trials_error ~ "0"))



#### Step 4: Make figure output of second smaller step ####
collapsed <- c()
for(i in 1:100) {
  # Start looking at the time series and plooting what they look like
  
  tree_name <- paste0(output_dir, "out_", i, ".trees")
  if(file.exists(tree_name)) { 
    
    tree <- read_annotated(filename = tree_name, format = "newick") 
    tree_plot <- ggtree(tree) + theme_tree2() + scale_x_continuous(breaks = seq(0, 20, 1))
    
    ts <- read.delim(paste0(output_dir, paste0("out_", i, ".timeseries")))
    
    if(subset(ts, date >5)$totalCases[1] == 0)
    {
      collapsed <- c(collapsed, i)
      next 
    } else { 
      
      plot_I <- ts %>% 
        pivot_longer(names_to = "case_type", values_to = "value", cols = c("totalI", "totalCases")) %>% 
        ggplot(aes(x = date, y = value, color = case_type)) + 
        geom_line() + theme_bw() + 
        labs(x = "date", y = "Value", color = "Case Type") +
        scale_color_manual(values = c("black", "red"), labels = c("Incidence", "Prevalence"))+
        scale_x_continuous(breaks = seq(1:20))
      
      ts$ar <- ts$totalI/ts$totalS*100
      
      
      tmrca <- read.csv(paste0(output_dir, paste0("out_", i, ".mrcaSeries")))
      colnames(tmrca) <- c("time", "tmrca")
      
      # The TMRCA plot -- also helps with quality control 
      plot_tmrca <- tmrca %>% 
        ggplot(aes(x = time, y = tmrca)) + geom_line() + 
        labs(x = "time (years)", y = "The Most Recent  \nCommon Ancestor") +
        scale_x_continuous(breaks = seq(1:20))
      
      # Antigen Series
      antigen_ts <- read_csv(paste0(output_dir, "out_", i, ".antigenicSamples"), col_names = FALSE)
      num_voc <- ncol(antigen_ts) - 1
      
      colnames(antigen_ts)[1] <- "time"
      colnames(antigen_ts)[2:ncol(antigen_ts)] <- paste0("voc_", seq(1:num_voc))
      
      
      mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(num_voc)
      antigen_ts <- antigen_ts %>% 
        pivot_longer(names_to = "voc", values_to = "numbers_sampled", -time, names_prefix = "voc_") 
      
      antigen_ts$voc <- factor(antigen_ts$voc, levels = seq(1:num_voc))
      
      plot_voc <- antigen_ts %>% 
        ggplot(aes(x = time, y = numbers_sampled, fill = voc)) + geom_bar(stat = "identity") +
        theme_bw() + labs(y = "Number of VOC Sampled", x = "time (years)") +
        scale_fill_manual(values = mycolors) + 
        scale_x_continuous(breaks = seq(1:20))
      
      # Put proportion of samples onto timeseries 
      voc_prop <- antigen_ts %>% 
        group_by(time) %>% 
        mutate(rel_freq = numbers_sampled/sum(numbers_sampled)) %>% 
        select(time, voc, rel_freq) 
      
      voc_prop <- voc_prop %>% 
        left_join(ts[, c("date", "totalI", "totalCases")], by = c("time" = "date"))  %>% 
        mutate(voc_I = rel_freq * totalCases) 
      
      # population from yaml files
      pop <- as.numeric(small_parameter_df[i, ]$initialNs)
      plot_voc_I <- voc_prop %>% 
        #filter(time > 10 & time < 19) %>% 
        ggplot(aes(x = time, y = voc_I/pop * 100000, fill = voc)) + 
        geom_area() + theme_bw() + scale_x_continuous(breaks = seq(1:20)) + 
        labs(y = "Incidence per 100,000", x = "Years") +
        scale_fill_manual(values = mycolors) + guides(fill = "none") +
        theme(text = element_text(size = 16)) 
      
      
      plot_turnover <- voc_prop %>% 
        #filter(time > 10 & time < 19) %>% 
        ggplot(aes(x = time, y = rel_freq, fill = voc)) + geom_area() + 
        scale_fill_viridis_d() + theme_bw() +
        labs(x = "Years", y = "Relative Frequency") + 
        scale_x_continuous(breaks = seq(1:20)) +
        scale_fill_manual(values = mycolors) + 
        theme(legend.position = "bottom") + 
        theme(legend.direction = "horizontal")
      
      plot_legend <- get_legend_35(plot_turnover) # + guides(fill=guide_legend(nrow=4)))
      
      # Put all together
      plot_col <- plot_grid(tree_plot, plot_tmrca, plot_I, 
                            plot_voc + theme(legend.position = "none"), 
                            plot_voc_I + theme(legend.position = "none"), 
                            plot_turnover + theme(legend.position = "none"), ncol = 2, align = "h", 
                            labels = "AUTO", label_x = .95, label_y = .95)
      
      plot_col <- plot_grid(plot_col, plot_legend, ncol = 1, rel_heights = c(.9, .1))  
      
      # make a plot grid consisting of two panels
      # now add the title
      title <- ggdraw() + 
        draw_label(
          paste0("Trial ", i, ""),
          fontface = 'bold',
          x = 0,
          hjust = 0
        ) +
        theme(
          # add margin on the left of the drawing canvas,
          # so title is aligned with left edge of first plot
          plot.margin = margin(0, 0, 0, 7)
        )
      final_plot <- plot_grid(
        title, plot_col,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.1, 1)
      )
      save_plot(filename = paste0(figure_dir, "trial_", i, ".png"), plot = final_plot, base_height = 8)
    }
  }
  else {
    print(paste0("Skipping", i, " - File not found"))
    next
  }
}




#### Step 5: Calculate evo/epi outcomes for the smaller set ####
small_outcomes <- data.frame()
for(i in trials_success) {
  # Start looking at the time series and plooting what they look like
  
  ts_name <- paste0(output_dir, "out_", i, ".timeseries")
    
  ts <- read.delim(paste0(output_dir, paste0("out_", i, ".timeseries")))
  ts$ar <- ts$totalI/ts$totalS*100
  average_ar <- mean(ts$ar)
  sd_ar <- sd(ts$ar)
  
  antigen_ts <- read_csv(paste0(output_dir, "out_", i, ".antigenicSamples"), col_names = FALSE)
  num_voc <- ncol(antigen_ts)-1
  turnover <- as.numeric(ifelse(which.max(antigen_ts[nrow(antigen_ts), 2:ncol(antigen_ts)]) == 1, 0, 1))
  
  # Number of peaks
  # Need to smooth 
  totalI_smooth <- loess(ts$totalI ~ ts$date, span = .05)
  ts$totalI_smooth <- predict(totalI_smooth)
  
  peaks <- findpeaks(ts$totalI_smooth, minpeakdistance = 15, nups = 5,  minpeakheight = 300, npeaks = 20)
  peak_heights <- peaks[,1]
  
  # code for checking heights
  peak_times <- data.frame(date = ts$date[peaks[, 2]])
  ts %>% 
    ggplot(aes(x = date, y = total_smooth)) + geom_line() + 
    geom_line(aes(x = date, y = totalI), color = "blue") +
    geom_vline(data = peak_times, aes(xintercept = date))
  
  npeaks <- nrow(peaks)
  sd_heights <- sd(peak_heights)    
  
  if(is.null(npeaks)) {
    npeaks <- 0
    sd_heights <- 0
  }
  
  
  run_profile <- data.frame(ar= average_ar, num_voc = num_voc, sd_epi_heights = sd_heights, num_epi_peaks = npeaks, turnover = turnover)
  small_outcomes <- rbind(small_outcomes, run_profile)
}

small_outcomes$trial <- trials_success
small_outcomes <- small_outcomes %>% left_join(small_parameter_df, by = c("trial" = "trial_num"))
small_outcomes$disease <- "synthetic"
table(small_outcomes$turnover)



#### Step 6: Now repeat for previous experiments ####
yaml_dir <- "../mutantigen_parallel/input_files/"

# Get a list of all YAML files in the directory
yaml_files <- list.files(yaml_dir, pattern = "\\.yml$", full.names = TRUE)

# Initialize an empty list to store data frames
parameter_list <- list()

# Loop through each YAML file
for (yaml_file in yaml_files) {
  # Load the YAML content into a list
  yaml_content <- yaml.load_file(yaml_file)
  
  # Convert the YAML content to a flat named vector (parameter names and values)
  flat_content <- unlist(yaml_content, recursive = TRUE)
  
  # Convert the vector to a data frame with parameter names as column headers
  parameter_df <- as.data.frame(t(flat_content), stringsAsFactors = FALSE)
  
  # Add a column for the file name
  parameter_df$file <- basename(yaml_file)
  
  # Append the data frame to the list
  parameter_list[[length(parameter_list) + 1]] <- parameter_df
}

large_parameter_df <- bind_rows(parameter_list)
large_parameter_df <- large_parameter_df %>% mutate(trial_num = as.numeric(str_extract(file,"[0-9]{1,5}")))

output_dir <- "../mutantigen_parallel_output/data/"
figure_dir <- "../mutantigen_parallel_output/figures/"

large_ran <- unique(as.numeric(str_extract(list.files(output_dir), pattern = "[0-9]{1,4}")))
large_success <- sort(as.numeric(str_extract(list.files(figure_dir), pattern = "[0-9]{1,4}")))

large_parameter_df <- large_parameter_df %>% 
  filter(trial_num %in% large_ran) %>% 
  mutate(category = case_when(trial_num %in% large_success ~ "1",
                              TRUE ~ "0"))

#### Step 7: Get larger set outcomes ####

large_outcomes <- data.frame()
for(i in large_success) {
  # Start looking at the time series and plooting what they look like
  
  ts_name <- paste0(output_dir, "out_", i, ".timeseries")
  
  ts <- read.delim(paste0(output_dir, paste0("out_", i, ".timeseries")))
  ts$ar <- ts$totalI/ts$totalS*100
  average_ar <- mean(ts$ar)
  sd_ar <- sd(ts$ar)
  
  antigen_ts <- read_csv(paste0(output_dir, "out_", i, ".antigenicSamples"), col_names = FALSE)
  num_voc <- ncol(antigen_ts)-1
  turnover <- as.numeric(ifelse(which.max(antigen_ts[nrow(antigen_ts), 2:ncol(antigen_ts)]) == 1, 0, 1))
  
  # Number of peaks
  # Need to smooth 
  totalI_smooth <- loess(ts$totalI ~ ts$date, span = .05)
  ts$totalI_smooth <- predict(totalI_smooth)
  
  peaks <- findpeaks(ts$totalI_smooth, minpeakdistance = 15, nups = 5,  minpeakheight = 300, npeaks = 20)
  peak_heights <- peaks[,1]
  
  # code for checking heights
  peak_times <- data.frame(date = ts$date[peaks[, 2]])
  ts %>% 
    ggplot(aes(x = date, y = total_smooth)) + geom_line() + 
    geom_line(aes(x = date, y = totalI), color = "blue") +
    geom_vline(data = peak_times, aes(xintercept = date))
  
  npeaks <- nrow(peaks)
  sd_heights <- sd(peak_heights)    
  
  if(is.null(npeaks)) {
    npeaks <- 0
    sd_heights <- 0
  }
  
  
  run_profile <- data.frame(ar= average_ar, num_voc = num_voc, sd_epi_heights = sd_heights, num_epi_peaks = npeaks, turnover = turnover)
  large_outcomes <- rbind(large_outcomes, run_profile)
}

large_outcomes$trial <- large_success
large_outcomes <- large_outcomes %>% left_join(large_parameter_df, by = c("trial" = "trial_num"))
large_outcomes$disease <- "synthetic"




###### VISUALIZATIONS & MODELS #####

# Question 1: What parameters ran and what parameters didn't
head(large_parameter_df)
head(small_parameter_df)
head(flu_parameter)

#### Step 8: Combine all parameter data frames into a single data frame ####

large_parameter_sub <- large_parameter_df[, c(which(colnames(large_parameter_df) %in% c("category", "initialI", names(param_ranges))))]
#large_parameter_sub <- large_parameter_sub[,-which(colnames(large_parameter_sub) == "epsilon_mut")]

small_parameter_sub <- small_parameter_df[, c(which(colnames(small_parameter_df) %in% c("category", "initialI", names(param_ranges))))]
#small_parameter_sub <- small_parameter_sub[,-which(colnames(small_parameter_sub) == "epsilon_mut")]

#flu_parameter_sub <- flu_parameter[, c(which(colnames(flu_parameter) %in% c("category", "epsilon", names(param_ranges))))]

super_parameter <- rbind(large_parameter_sub, small_parameter_sub)
#super_parameter <- rbind(small_parameter_sub, flu_parameter_sub)

cols_to_convert <- names(super_parameter)[sapply(super_parameter, is.character)] # Identify character columns

super_parameter[cols_to_convert] <- lapply(super_parameter[cols_to_convert], as.numeric)
super_parameter$initialI_prop   <- super_parameter$initialI/super_parameter$initialNs
super_parameter$category <- factor(super_parameter$category)
colnames(super_parameter)
super_parameter <- super_parameter[,-which(colnames(super_parameter) == "initialI")]


# Visualization of what 
ggpairs(super_parameter,                # Data frame
       columns = c(1:9, 11),        # Columns
       aes(color = as.factor(category),  # Color by group (cat. variable)
           alpha = 0.5))  




#### Step 9: Put together for model ####
library(caret)
set.seed(022924)

# split up data into 80% training and 20% testing
levels(super_parameter$category) <- make.names(levels(super_parameter$category))  # Ensures valid names
print(levels(super_parameter$category))  # Check new levels
super_parameter$category <- relevel(super_parameter$category, ref = "X0")


training_ids <- sample(x = seq(from = 1, to = nrow(super_parameter)), size = .80 * nrow(super_parameter), replace = FALSE)

train_df <- super_parameter[training_ids, ]
test_df <- super_parameter[-training_ids,]

 sensitivity_metric <- function(data, lev = NULL, model = NULL) {
   cm <- confusionMatrix(data$pred, data$obs, positive = "X1")
   #specificity <- cm$byClass["Specificity"]
   #c(Specificity = specificity)
   sensitivity <- cm$byClass["Sensitivity"]
   c(Sensitivity = sensitivity)
 }


F1_metric <- function(data, lev = NULL, model = NULL) {
  if (!all(c("pred", "obs") %in% colnames(data))) {
    stop("Missing 'pred' or 'obs' column in the data")
  }
  
  precision <- posPredValue(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  recall <- sensitivity(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  
  # Ensure F1 is not NA
  if (is.na(precision) | is.na(recall) | (precision + recall) == 0) {
    F1 <- 0
  } else {
    F1 <- 2 * (precision * recall) / (precision + recall)
  }
  
  return(c(F1 = F1))
}

F1_weighted_metric <- function(data, lev = NULL, model = NULL, beta = 2) {
  if (!all(c("pred", "obs") %in% colnames(data))) {
    stop("Missing 'pred' or 'obs' column in the data")
  }
  
  precision <- posPredValue(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  recall <- sensitivity(data$pred, data$obs, positive = "X1", na.rm = TRUE)
  
  # Ensure F1 is not NA
  if (is.na(precision) | is.na(recall) | (beta^2 * precision + recall) == 0) {
    F_beta <- 0
  } else {
    F_beta <- (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall)
  }
  
  return(c(F1_beta = F_beta))
}


# Define trainControl for 5-fold cross-validation

###### GBM ####
train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = sensitivity_metric,
  sampling = "smote"
)

trainGrid <- expand.grid(n.trees = c(15, 20, 25),
                         interaction.depth = c(2, 3, 5), 
                         shrinkage = c(0.1),
                         n.minobsinnode = 20)

gbm_model_sensitivity <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
#  metric = "F1",
)

gbm_model_sensitivity
varImp(gbm_model, scale = TRUE)

# Evaluate training data model performance
gbm_train_predictions <- predict(gbm_model_sensitivity, train_df)
gbm_train_cm <- confusionMatrix(gbm_train_predictions, train_df$category, positive = "X1")

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_sensitivity, test_df)
gbm_test_cm_sensitivity <- confusionMatrix(gbm_test_predictions, test_df$category, positive = "X1")
#saveRDS(gbm_model_sensitivity, file = "../mutantigen_experiment/gbm_sensitivity.RDS")



# Look at where it gets answers wrong
test_df$predictions <- gbm_test_predictions
test_df$result <- test_df$predictions == test_df$category

# Visualization of what 
ggpairs(test_df,                # Data frame
        columns = c(1:9, 11:12),        # Columns
        aes(color = as.factor(result),  # Color by group (cat. variable)
            alpha = 0.5))  

# Looks like we are getting lambda predictions wrong -- in which direction?
actual_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = category)) + geom_density(alpha = .5)

predicted_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = predictions)) + geom_density(alpha = .5)

accuracy_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = result)) + geom_density(alpha = .5)


# when lambda is low, there is a tendency to think it will finish, but it may not always
# Not missing anything
plot_grid(actual_split, predicted_split, accuracy_split, ncol = 1)



##### Testing for other control measures ####

train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = F1_metric,
  sampling = "smote"
)

trainGrid <- expand.grid(n.trees = c(20, 35, 50, 100, 500),
                         interaction.depth = c(3, 5, 10, 20), 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)
gbm_model_f1 <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
  #  metric = "F1",
)
saveRDS(gbm_model_f1, file = "../mutantigen_experiment/gbm_f1_step1.RDS")




gbm_train_predictions <- predict(gbm_model_f1, train_df)
gbm_train_cm <- confusionMatrix(gbm_train_predictions, train_df$category, positive = "X1")

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_f1, test_df)
gbm_test_cm_f1 <- confusionMatrix(gbm_test_predictions, test_df$category, positive = "X1")
print(gbm_test_cm)

test_df$predictions <- gbm_test_predictions
test_df$result <- test_df$predictions == test_df$category

# Visualization of what 
ggpairs(test_df,                # Data frame
        columns = c(1:9, 11:12),        # Columns
        aes(color = as.factor(result),  # Color by group (cat. variable)
            alpha = 0.5))  

# Looks like we are getting lambda predictions wrong -- in which direction?
actual_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = category)) + geom_density(alpha = .5)

predicted_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = predictions)) + geom_density(alpha = .5)

accuracy_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = result)) + geom_density(alpha = .5)


# when lambda is low, there is a tendency to think it will finish, but it may not always
# Not missing anything
plot_grid(actual_split, predicted_split, accuracy_split, ncol = 1)



train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = function(data, lev, model) F1_weighted_metric(data, lev, model, beta = 2),  # More weight on recall
  sampling = "smote"
)

trainGrid <- expand.grid(n.trees = c(20, 35, 50, 100, 500),
                         interaction.depth = c(3, 5, 10, 20), 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)

gbm_model_f1_weighted <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
    metric = "F1_beta"
)


gbm_train_predictions <- predict(gbm_model_f1_weighted, train_df)
gbm_train_cm <- confusionMatrix(gbm_train_predictions, train_df$category, positive = "X1")

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_f1_weighted, test_df)
gbm_test_cm_weighted <- confusionMatrix(gbm_test_predictions, test_df$category, positive = "X1")

## Looking at my three different metrics
results <- data.frame(sensitivity = gbm_test_cm_sensitivity$byClass,
                      f1 = gbm_test_cm_f1$byClass,
                      fbeta_2 = gbm_test_cm_weighted$byClass)
# going with F1, seems to be a good balance





### Step 10: Now visualize model if there is turnover or not ####
large_outcomes_sub <- large_outcomes[, c(which(colnames(large_outcomes) %in% c("turnover", "initialI", names(param_ranges))))]
#large_outcomes_sub <- large_outcomes_sub[,-which(colnames(large_outcomes_sub) == "epsilon_mut")]

small_outcomes_sub <- small_outcomes[, c(which(colnames(small_outcomes) %in% c("turnover", "initialI", names(param_ranges))))]
#small_outcomes_sub <- small_outcomes_sub[,-which(colnames(small_outcomes_sub) == "epsilon_mut")]

#flu_outcomes$turnover <- "1"
#flu_outcomes_sub <- flu_outcomes[, c(which(colnames(flu_outcomes) %in% c("turnover", "epsilon", names(param_ranges))))]

#super_outcomes <- rbind(large_outcomes_sub, small_outcomes_sub, flu_outcomes_sub[1,])
super_outcomes <- rbind(large_outcomes_sub, small_outcomes_sub)

#super_outcome <- rbind(small_outcome_sub, flu_outcome_sub)

cols_to_convert <- names(super_outcomes)[sapply(super_outcomes, is.character)] # Identify character columns
super_outcomes[cols_to_convert] <- lapply(super_outcomes[cols_to_convert], as.numeric)
super_outcomes$turnover <- factor(super_outcomes$turnover)
super_outcomes$initialI_prop   <- super_outcomes$initialI/super_outcomes$initialNs
super_outcomes <- super_outcomes[,-which(colnames(super_outcomes) == "initialI")]

#super_outcomes$epsilon_pop <- super_outcomes$epsilon/super_outcomes$initialNs
#super_outcomes$r0 <- super_outcomes$beta/super_outcomes$nu

# Visualization of what 
ggpairs(super_outcomes,                # Data frame
        aes(color = as.factor(turnover),  # Color by group (cat. variable)
            alpha = 0.5))  

#### Step 11: Put together for model ####
training_ids <- sample(x = seq(from = 1, to = nrow(super_outcomes)), size = .80 * nrow(super_outcomes), replace = FALSE)

levels(super_outcomes$turnover) <- make.names(levels(super_outcomes$turnover))  # Ensures valid names
print(levels(super_outcomes$turnover))  # Check new levels

super_outcomes$turnover <- relevel(super_outcomes$turnover, ref = "X0")

train_df <- super_outcomes[training_ids, ]
test_df <- super_outcomes[-training_ids,]

summary(train_df)
summary(test_df)


train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = sensitivity_metric,
  sampling = "smote"
)

trainGrid <- expand.grid(n.trees = c(15, 20, 25),
                         interaction.depth = c(2, 3, 5), 
                         shrinkage = c(0.1),
                         n.minobsinnode = c(5, 10, 20))

gbm_model_sensitivity <- train(
  turnover ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
  #  metric = "F1",
)
gbm_model_sensitivity
varImp(gbm_model_sensitivity)

# Evaluate training data model performance
predictions <- predict(gbm_model_sensitivity, train_df)
confusion_matrix <- confusionMatrix(predictions, train_df$turnover, positive = "X1")


# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_sensitivity, test_df)
# Specificity is low straight out of the box ~ 20%
gbm_test_sensitivity <- confusionMatrix(gbm_test_predictions, test_df$turnover, positive = "X1")


# Look at where it gets answers wrong
test_df$predictions <- gbm_test_predictions
test_df$result <- test_df$predictions == test_df$turnover

# Visualization of what 
ggpairs(test_df,                # Data frame
        columns = c(2:12),        # Columns
        aes(color = as.factor(result),  # Color by group (cat. variable)
            alpha = 0.5))  


# Looks like we are getting lambda predictions wrong -- in which direction?
actual_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = category)) + geom_density(alpha = .5)

predicted_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = predictions)) + geom_density(alpha = .5)

accuracy_split <- test_df %>% 
  ggplot(aes(x = lambda, fill = result)) + geom_density(alpha = .5)


# when lambda is low, there is a tendency to think it will finish, but it may not always
# Not missing anything
plot_grid(actual_split, predicted_split, accuracy_split, ncol = 1)






#### Other Metrics ####

train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = F1_metric,
  sampling = "smote"
)

trainGrid <- expand.grid(n.trees = c(20, 35, 50, 100, 500),
                         interaction.depth = c(3, 5, 10, 20), 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)
gbm_model_f1 <- train(
  turnover ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
  #  metric = "F1",
)


gbm_train_predictions <- predict(gbm_model_f1, train_df)
gbm_train_cm <- confusionMatrix(gbm_train_predictions, train_df$turnover, positive = "X1")

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_f1, test_df)
gbm_test_cm_f1 <- confusionMatrix(gbm_test_predictions, test_df$turnover, positive = "X1")



# Look at where it gets answers wrong
test_df$predictions <- gbm_test_predictions
test_df$result <- test_df$predictions == test_df$turnover

# Visualization of what 
ggpairs(test_df,                # Data frame
        columns = c(2:12),        # Columns
        aes(color = as.factor(result),  # Color by group (cat. variable)
            alpha = 0.5))  


# lambda and mean antigenic size
#saveRDS(gbm_model_f1, "../mutantigen_experiment/gbm_f1_step2.RDS")





train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  classProbs = TRUE,
  verboseIter = TRUE,  # Show progress
  summaryFunction = function(data, lev, model) F1_weighted_metric(data, lev, model, beta = 2),  # More weight on recall
  sampling = "smote"
)

trainGrid <- expand.grid(n.trees = c(20, 35, 50, 100, 500),
                         interaction.depth = c(3, 5, 10, 20), 
                         shrinkage = 0.1,
                         n.minobsinnode = 20)

gbm_model_f1_weighted <- train(
  turnover ~ .,         # Formula (target ~ predictors)
  data = train_df,          # Data frame
  method = "gbm",      # Logistic Regression
  trControl = train_control,
  tuneGrid = trainGrid,
  #  weights = ifelse(train_df$category == "1", 2, 1),
  tuneLength = 5,
  metric = "F1_beta"
)


gbm_train_predictions <- predict(gbm_model_f1_weighted, train_df)
gbm_train_cm <- confusionMatrix(gbm_train_predictions, train_df$turnover, positive = "X1")

# Evaluate testing data model performance
gbm_test_predictions <- predict(gbm_model_f1_weighted, test_df)
gbm_test_cm_weighted <- confusionMatrix(gbm_test_predictions, test_df$turnover, positive = "X1")

## Looking at my three different metrics
smote_sampling_results <- data.frame(sensitivity = gbm_test_sensitivity$byClass,
                      f1 = gbm_test_cm_f1$byClass,
                      fbeta_2 = gbm_test_cm_weighted$byClass)
saveRDS(gbm_model_f1_weighted, "../mutantigen_experiment/gbm_fbeta_step2.RDS")

# training on sensitivity with down sampling 




######### GRAVEYARD ##################




# Put canonical flu on here 

run_outcome_sub <- run_profile_df[, c(which(colnames(run_profile_df) %in% c("turnover", "epsilon", names(param_ranges))))]
cols_to_convert <- names(run_outcome_sub)[sapply(run_outcome_sub, is.character)] # Identify character columns

run_outcome_sub[cols_to_convert] <- lapply(run_outcome_sub[cols_to_convert], as.numeric)
run_outcome_sub$turnover <- factor(run_outcome_sub$turnover)
run_outcome_sub$epsilon_pop <- run_outcome_sub$epsilon/run_outcome_sub$initialNs
run_outcome_sub$r0 <- run_outcome_sub$beta/run_outcome_sub$nu

# This figure shows impacts of differnt evo parameters on turnover
# Also shows where we have some "non" runs

# Most apparent differences are in beta, lambdaantigenic, potential epsilon_pop, lambd

ggpairs(run_outcome_sub,       # Data frame
        columns = c(2:7, 9:10, 12:13),        # Columns
        aes(color = as.factor(turnover),  # Color by group (cat. variable)
            alpha = 0.5))      # Transparency

ggpairs(run_profile_df,
        columns = 1:4,
        aes(color = factor(turnover),
            alpha = 0.5))


# run_profile_df %>% 
#   select(num_epi_peaks, turnover, colnames(params_run)) %>% 
#   pivot_longer(names_to = "outcome", values_to = "value", -c(turnover, num_epi_peaks)) %>% 
#   mutate(value = as.numeric(value)) %>% 
#   ggplot(aes(x = value, y = num_epi_peaks, color = as.factor(turnover))) + 
#   geom_point() + facet_wrap(~outcome, scales = "free")


run_profile_df %>% 
  mutate(disease = "sample") %>% 
  select(num_epi_peaks, turnover, colnames(params_run)) %>% 
  pivot_longer(names_to = "outcome", values_to = "value", -c(turnover, num_epi_peaks)) %>% 
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(x = value, y = num_epi_peaks, color = as.factor(turnover))) + 
  geom_point() + facet_wrap(~outcome, scales = "free")


flu_profile_df$turnover <- factor(flu_profile_df$turnover)
run_profile_df$turnover <- factor(run_profile_df$turnover)

# showing difference of flu vs disease 
run_profile_df %>% 
  mutate(disease = "sample") %>% 
  bind_rows(flu_profile_df) %>% 
  select(ar, num_voc, sd_epi_heights, num_epi_peaks, disease)  %>% 
  pivot_longer(names_to = "outcome", values_to = "value", -disease) %>% 
  ggplot(aes(x = value, color = disease)) + geom_density() + facet_wrap(~outcome, scales = "free")



combined <- run_profile_df %>% 
  mutate(disease = "sample") %>% 
  bind_rows(flu_profile_df)

combined_sub <- combined[, c(which(colnames(combined) %in% c("disease", "epsilon", names(param_ranges))))]
cols_to_convert <- names(combined_sub)[sapply(combined_sub, is.character)] # Identify character columns
combined_sub[cols_to_convert[-length(cols_to_convert)]] <- lapply(combined_sub[cols_to_convert[-length(cols_to_convert)]], as.numeric)
combined_sub <- combined_sub[, -which(colnames(combined_sub) == "epsilon_mut")]
combined_sub$epsilon_pop = combined_sub$epsilon/combined_sub$initialNs
plot_parameter <- ggpairs(combined_sub,                    # Data frame
        columns = c(1:6,8:9, 11),        # Columns
        aes(color = disease,  # Color by group (cat. variable)
            alpha = 0.5))      # Transparency

cowplot::save_plot(plot = plot_parameter, 
                   filename = "../mutantigen_parallel_output/figures/parameter_distributions.png", 
                   base_height = 8)



combined_sub    <- combined[, c(which(colnames(combined) %in% c("disease", "ar", "num_voc", "sd_epi_heights", "num_epi_peaks")))]
cols_to_convert <- names(combined_sub)[sapply(combined_sub, is.character)] # Identify character columns
combined_sub[cols_to_convert[-length(cols_to_convert)]] <- lapply(combined_sub[cols_to_convert[-length(cols_to_convert)]],
                                                                  as.numeric) 
plot_outcome <- ggpairs(combined_sub,                    # Data frame
                          columns = c(1:4),        # Columns
                          aes(color = disease,  # Color by group (cat. variable)
                              alpha = 0.5))      # Transparency

cowplot::save_plot(plot = plot_outcome, 
                   filename = "../mutantigen_parallel_output/figures/outcome_distributions.png", 
                   base_height = 8)



### Look at the parameters run compared to what is missing: 

trials_success <- sort(as.numeric(str_extract(list.files(figure_dir), pattern = "[0-9]{1,3}")))
trials_error   <- seq(1:100)[-trials_success]

# What was filtered out or had an error
# What didn't run/
final_df <- final_df %>% 
  mutate(trialNum = as.numeric(str_extract(file, pattern = "[0-9]{1,5}"))) %>% 
  mutate(category = case_when(trialNum %in% trials_success ~ "1",
                              trialNum %in% trials_error ~ "0"))



final_df %>% 
  select(trialNum, initialNs, demeAmplitudes, lambdaAntigenic, meanAntigenicSize, 
         lambda, mutCost, beta, nu, epsilon_mut, initialI, category) %>% 
  pivot_longer(names_to = "param", values_to = "value", cols = -c(trialNum, category))  %>% 
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(x = value, color = category)) + 
  geom_density() +
  #geom_histogram(position = "dodge", stat = "count") + 
  facet_wrap(~param,scales = "free")


library(caret)
set.seed(022924)

# Define trainControl for 5-fold cross-validation
train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 3,         # Number of folds
  verboseIter = TRUE  # Show progress
)

df <- final_df[, c(which(colnames(final_df) %in% c("epsilon", "category", names(param_ranges))))]
cols_to_convert <- names(df)[sapply(df, is.character)] # Identify character columns
df[cols_to_convert] <- lapply(df[cols_to_convert], as.numeric)

df$category <- as.factor(df$category)
model <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = df,          # Data frame
  method = "rf",      # Logistic Regression
  trControl = train_control,
  ntree = 40
)

# Print model results
print(model)


# Make predictions on the same data 
predictions <- predict(model, df)

# Evaluate model performance
confusion_matrix <- confusionMatrix(predictions, df$category)
print(confusion_matrix)

svm_model <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = df,          # Data frame
  method = "svmRadial",      # Logistic Regression
  preProcess = c("center", "scale"),
  trControl = train_control,
  tuneLength = 5
)

predictions <- predict(svm_model)

# Evaluate model performance
confusion_matrix <- confusionMatrix(predictions, df$category)
print(confusion_matrix)


### Now test each of these models are the original runs
yaml_dir <- "../mutantigen_parallel/input_files/"

# Get a list of all YAML files in the directory
yaml_files <- list.files(yaml_dir, pattern = "\\.yml$", full.names = TRUE)

# Initialize an empty list to store data frames
parameter_list <- list()

# Loop through each YAML file
for (yaml_file in yaml_files) {
  # Load the YAML content into a list
  yaml_content <- yaml.load_file(yaml_file)
  
  # Convert the YAML content to a flat named vector (parameter names and values)
  flat_content <- unlist(yaml_content, recursive = TRUE)
  
  # Convert the vector to a data frame with parameter names as column headers
  parameter_df <- as.data.frame(t(flat_content), stringsAsFactors = FALSE)
  
  # Add a column for the file name
  parameter_df$file <- basename(yaml_file)
  
  # Append the data frame to the list
  parameter_list[[length(parameter_list) + 1]] <- parameter_df
}

# Combine all parameter data frames into a single data frame
testing_df <- bind_rows(parameter_list)
testing_df <- testing_df %>% mutate(trial_num = as.numeric(str_extract(file,"[0-9]{1,5}")))

output_dir <- "../mutantigen_parallel_output/data/"
figure_dir <- "../mutantigen_parallel_output/figures/"

testing_ran <- unique(as.numeric(str_extract(list.files(output_dir), pattern = "[0-9]{1,4}")))
length(testing_ran)
testing_success <- sort(as.numeric(str_extract(list.files(figure_dir), pattern = "[0-9]{1,4}")))
length(testing_success)

292/(1854+292)

testing_df <- testing_df %>% 
  filter(trial_num %in% testing_ran) %>% 
  arrange(trial_num)

testing_df <- testing_df %>% 
  mutate(category = case_when(trial_num %in% testing_success ~ "1",
                              trial_num %in% testing_ran ~ "0"))


test_df <- testing_df[, c(which(colnames(testing_df) %in% c("epsilon", "category", names(param_ranges))))]
cols_to_convert <- names(test_df)[sapply(test_df, is.character)] # Identify character columns
test_df[cols_to_convert] <- lapply(test_df[cols_to_convert], as.numeric)


test_predictions <- predict(model, newdata = test_df)

# Evaluate model performance
test_df$category <- as.factor(test_df$category)
confusion_matrix <- confusionMatrix(test_predictions, test_df$category)
print(confusion_matrix)


# Evaluate model performance
svm_test_predictions <- predict(svm_model, newdata = test_df)
test_df$category <- as.factor(test_df$category)
confusion_matrix <- confusionMatrix(test_predictions, test_df$category)
print(confusion_matrix)



# Model is overfiting (probably need more training data)
# I have aboout 2000 simulations -- save 200 of these for testing. 

# Pull in all the new simulations
# Pull in some of the flu simulations
# Pull in some of the original simulations
# Test on rest of the original simulations

#flu_sub <- flu_profile_df[sample(1:30, size = 5, replace = 5),]
flu_parameter$category <- 1
flu_parameter <- flu_parameter[1, c(which(colnames(flu_parameter) %in% c("epsilon", "category", names(param_ranges))))] 
#flu_sub <- flu_sub %>% select(-disease)

validation_nums <- sample(1:nrow(testing_df), size = 1700, replace = F)
validation_df   <- testing_df[validation_nums, c(which(colnames(testing_df) %in% c("epsilon", "category", names(param_ranges))))]
validation_df <- subset(validation_df, select = -epsilon_mut)

holdout_df  <- testing_df[-validation_nums,] 

final_df <-  final_df[, c(which(colnames(final_df) %in% c("epsilon", "category", names(param_ranges))))]
final_df <- subset(final_df, select = -epsilon_mut)

# combined set
super_df <- rbind(flu_parameter, validation_df, final_df)


cols_to_convert <- names(super_df)[sapply(super_df, is.character)] # Identify character columns
super_df[cols_to_convert] <- lapply(super_df[cols_to_convert], as.numeric)

super_df$category <- as.factor(super_df$category)


set.seed(022924)

custom_metric <- function(data, lev = NULL, model = NULL) {
  cm <- confusionMatrix(data$pred, data$obs, positive = "1")
  #specificity <- cm$byClass["Specificity"]
  #c(Specificity = specificity)
  sensitivity <- cm$byClass["Sensitivity"]
  c(Sensitivity = sensitivity)
}
# Define trainControl for 5-fold cross-validation

train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  verboseIter = TRUE,  # Show progress
  summaryFunction = custom_metric
)



model <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = super_df,          # Data frame
  method = "rf",      # Logistic Regression
  trControl = train_control,
  weights = ifelse(super_df$category == "positive_class", 2,1),
  metric = "Sensitivity"
)

svm_model <- train(
  category ~ .,         # Formula (target ~ predictors)
  data = super_df,          # Data frame
  method = "svmRadial",      # Logistic Regression
  preProcess = c("center", "scale"),
  trControl = train_control,
  tuneLength = 5
)


# Print model results
print(model)

# Make predictions on the same data 
predictions <- predict(model, super_df)

# Evaluate model performance
confusion_matrix <- confusionMatrix(predictions, super_df$category)
print(confusion_matrix)


# predict on hold out
cols_to_convert <- names(holdout_df)[sapply(holdout_df, is.character)] # Identify character columns
holdout_df[cols_to_convert] <- lapply(holdout_df[cols_to_convert], as.numeric)
holdout_df$category <- as.factor(holdout_df$category)

holdout_predictions <- predict(model, newdata = holdout_df)
confusionMatrix(holdout_predictions, holdout_df$category)


holdout_predictions <- predict(svm_model, newdata = holdout_df)
confusionMatrix(holdout_predictions, holdout_df$category)










#### For the initial runs I have -- 
# can we start predicting which ones are going to turnover (be good) and which ones aren't?

df <- run_outcome_sub[, c(which(colnames(run_outcome_sub) %in% c("epsilon", "turnover", names(param_ranges))))]
df$turnover <- as.factor(df$turnover)

library(caret)
set.seed(022924)

# Define trainControl for 5-fold cross-validation
train_control <- trainControl(
  method = "cv",      # Cross-validation
  number = 5,         # Number of folds
  verboseIter = TRUE  # Show progress
)

# Train the classifier (e.g., GLM)
model <- train(
  turnover ~ .,         # Formula (target ~ predictors)
  data = df,          # Data frame
  method = "glm",      # Logistic Regression
  family = "binomial",
  trControl = train_control
)

# Print model results
print(model)

# Make predictions on the same data (for demonstration purposes)
predictions <- predict(model, df)

# Evaluate model performance
predictions <- predict(model, df)
confusion_matrix <- confusionMatrix(predictions, df$turnover)
print(confusion_matrix)



# Train the classifier (e.g., randomForest) - clear winner 
df$r0 = df$beta/df$nu

# Can we remove any for redundancy?
cor(df[, 2:ncol(df)])
which(colnames(df) == "epsilon")

model <- train(
  num_voc ~ .,         # Formula (target ~ predictors)
  data = df[, -8],          # Data frame
  method = "rf",      # Logistic Regression
  trControl = train_control,
)

# Print model results
print(model)


# Make predictions on the same data 
predictions <- predict(model, df)

# Evaluate model performance
predictions <- predict(model, df)
confusion_matrix <- confusionMatrix(predictions, df$turnover)
print(confusion_matrix)

importance <- varImp(model)
print(importance)


# Repeat for some of the other parameters, like attack rates
df <- run_profile_df[, c(which(colnames(run_profile_df) %in% 
                           c("num_voc", "epsilon", names(param_ranges))))]
cols_to_convert <- names(df)[sapply(df, is.character)] # Identify character columns
df[cols_to_convert] <- lapply(df[cols_to_convert], as.numeric)
                                                                  
df$r0 <- df$beta/df$nu


# Can we remove any for redundancy?
cor(df[, 2:ncol(df)])
which(colnames(df) =="epsilon")


str(df$num_voc)

# Define trainControl for 5-fold cross-validation
train_control <- trainControl(method = "cv", number = 5, verboseIter = TRUE)

param_grid <- expand.grid(mtry = c(2, 3, 4), ntree = c(10, 20, 30))
param_grid <- data.frame(param_grid)

tunegrid <- expand.grid(.mtry = c(2,3,4))
modellist <- list()

#train with different ntree parameters
for (max in c(20, 40,80, 160)){
  set.seed(123)
  fit <- train(num_voc ~., 
               data = df[,-8],
               method = 'rf',
               tuneGrid = tunegrid,
               trControl = train_control,
               max
               ntree = ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}

#Compare results
results <- resamples(modellist)
summary(results)



model <- caret::train(
  num_voc ~ .,           # Formula (target ~ predictors)
  data = df[, -8],       # Data frame
  method = "rf",       # Random Forest
  tuneGrid = tunegrid,
  trControl = train_control
)

# Print model results
print(model)


# Make predictions on the same data 
predictions <- predict(model, df)


actual <- df$num_voc
mse <- mean((predictions - actual)^2)  # Mean Squared Error
cat("Mean Squared Error:", mse, "\n")

results <- data.frame(predictions = predictions, actual = actual)
results %>% 
  ggplot(aes(x = predictions, y = actual)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) + theme_bw()


summary(model)

importance <- varImp(model)
print(importance)








# Train the classifier (e.g., svmRadial)
model <- train(
  num_voc ~ .,         # Formula (target ~ predictors)
  data = df,          # Data frame
  method = "svmRadial",      # Logistic Regression
  preProcess = c("center", "scale"),
  trControl = train_control,
  tuneLength = 5
)

# Print model results
# Make predictions on the same data (for demonstration purposes)
predictions <- predict(model, df)

# Evaluate model performance
predictions <- predict(model, df)
confusion_matrix <- confusionMatrix(predictions, df$turnover)
print(confusion_matrix)













