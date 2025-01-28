# This script performs the following tasks:
# splits data into training and testing
# defines the structured and unstructured models
# plots the time series fits and projections
# calculates MAEs

# Code by Julie A. Spencer and Scottie Alexander
# ============================================================================ #
# Define the function to create cumulative sum, find the threshold, and split data
split_data <- function(df, threshold = 300) {
  
  # Create the cumulative sum column 'sum'
  cum_cases <- cumsum(df$all)
  # Find the first time step where the cumulative sum is >= threshold
  threshold_time <- df$time[min(which(cum_cases >= threshold))]
  # Create the training dataset (from the first time step to the threshold time step)
  training_data <- df[df$time <= threshold_time, ]
  # Remove the "sum" column to prep for fitting
  # training_data <- select(df, -sum)
  # Create the testing dataset (after the threshold time step)
  testing_data <- df[df$time > threshold_time, ]
  # Remove the "sum" column to prep for fitting
  #testing_data <- select(df, -sum)
  return (list(training_data=training_data, testing_data=testing_data))
}
# ============================================================================ #
# Function that evaluates the I1-I2 model
exponential_2state <- function (t, x, params) {
  # Extract the state variables
  I1 = x['I1'] # state one
  I2 = x['I2'] # state two
  N = I1 + I2
  # Extract the parameter values
  beta11 = params[1]
  beta12 = params[2]
  beta21 = params[3]
  beta22 = params[4]
  gamma = params[5]
  # Evaluate the ordinary differential equations at time t
  dI1dt = I1*beta11 +I2*beta21 - I1*gamma
  dI2dt = I2*beta22 + I1*beta12 - I2*gamma
  # Combine into a single vector
  dxdt = c(dI1dt, dI2dt)
  # Return as a list (required by ODE solver)
  return(list(dxdt))
}
# ============================================================================ #
# Function that evaluates the simple (unstructured) exponential model
exponential_1state <- function (t, x, params) {
  # Extract the state variable
  X = x['X'] # state one
  # Extract the parameter values
  r = params[1]
  X0 = params[2]
  # Evaluate the ordinary differential equations at time t
  dX = r*X
  # Combine into a single vector
  dxdt = c(dX)
  # Return as a list (required by ODE solver)
  return(list(dxdt))
}
# ============================================================================ #
# struct_data - a dataframe (not list thereof) of simulated structured data
plot_single_projection <- function(struct_data, struct_params, struct_model,
                                   unstruct_data, unstruct_params, unstruct_model) {
  
  struct_split <- split_data(struct_data)
  unstruct_split <- split_data(unstruct_data)
  
  # -------------------------------------------------------------------------- #
  tsfit <- fitode(struct_model, data = struct_split$training_data,
                     start = struct_params,
                     tcol = "time"
  )
  
  cf <- tsfit@coef
  
  # start and end time for model interpolation
  t0 <- struct_split$training_data$time[1]
  t1 <- tail(struct_split$training_data$time, n=1)

  # these "predictions" are really just interpolating w/in the training data
  # using the fitted model
  pred <- predict(tsfit, times=seq(t0, t1, length.out=100))
  solns <- data.frame(time=pred$adults$times, adults=pred$adults$estimate,
                      children=pred$children$estimate)
  
  # training / testing split time, and new end time that includes extrapolation
  # period (i.e. true prediction / projection)
  t_split <- struct_split$testing_data$time[1] - 1
  t1 <- tail(struct_split$testing_data$time, n=1)
  
  # projections at high temporal resolution for plotting (includes training and
  # testing periods)
  pred <- predict(tsfit, times=seq(t0, t1, length.out=100))
  projected_solns_plot <- data.frame(time=pred$adults$times,
                                     adults=pred$adults$estimate,
                                     children=pred$children$estimate)
  
  # projections for calculating MAE (~ same time resolution as the data)
  pred <- predict(tsfit, times=seq(struct_split$training_data$time[1], t1))
  projected_solns <- data.frame(time=pred$adults$times, adults=pred$adults$estimate,
                                children=pred$children$estimate)
  
  # true predictions / projections for calculating MAE should only include times
  # present w/in the testing data
  projected_solns <- projected_solns %>% filter(time %in% struct_split$testing_data$time)
  
  # -------------------------------------------------------------------------- #
  fig <- ggplot(struct_data, aes(x = time)) +
    geom_point(aes(y = adults, color = "adults-cases")) +  # Plot for data1
    geom_point(aes(y = children, color = "children-cases")) +  # Plot for data2
    labs(title = "", x = "Time", y = "Cases") +
    geom_line(data = solns, aes(x = time, y = adults, color = "adults-fitted")) +  # Fitted line for Data1
    geom_line(data = solns, aes(x = time, y = children, color = "children-fitted")) +  # Fitted line for Data2
    geom_line(data = projected_solns_plot, aes(x = time, y = adults, color = "adults-projected"),linetype="dashed") + # I1 forecast
    geom_line(data = projected_solns_plot, aes(x = time, y = children, color = "children-projected"),linetype = "dashed") + #I2 forecast
    scale_y_log10() +
    geom_vline(xintercept=t_split + 0.1, linetype="dotted", color="gray22", linewidth = 0.7) +
    scale_color_manual(name = "Legend",
                       values = c("adults-cases" = "blue",
                                  "children-cases" = "orange",
                                  "adults-fitted" = "blue2",
                                  "children-fitted" = "darkorange2",
                                  "adults-projected" = "blue2",
                                  "children-projected" = "darkorange2")) +
    theme_minimal()
  # -------------------------------------------------------------------------- #
  tsfit_un <- fitode(unstruct_model, data = unstruct_split$training_data,
                       start = unstruct_params,
                       tcol = "time"
  )
  
  # start and end time for model interpolation
  t0 <- unstruct_split$training_data$time[1]
  t1 <- tail(unstruct_split$training_data$time, n=1)

  solns_one <- predict(tsfit_un, times=seq(t0, t1, length.out=100))$all
  colnames(solns_one) <- c("time","X")
  
  # training / testing split time, and new end time that includes extrapolation
  # period (i.e. true prediction / projection)
  t_split <- unstruct_split$testing_data$time[1] - 1
  t1 <- tail(unstruct_split$testing_data$time, n=1)
  
  # projections at high temporal resolution for plotting
  projected_solns_one_plot <- predict(tsfit_un, times=seq(t0, t1, length.out=100))$all
  colnames(projected_solns_one_plot) <- c("time","X")
  
  # projections for calculating MAE (~ same time resolution as the data)
  pred <- predict(tsfit_un, times=seq(unstruct_split$training_data$time[1], t1))
  projected_solns_one <- data.frame(time=pred$all$times, X=pred$all$estimate)
  
  projected_solns_one <- projected_solns_one %>% filter(time %in% unstruct_split$testing_data$time)
  # -------------------------------------------------------------------------- #
  fig_un <- ggplot(unstruct_data, aes(x = time)) +
    geom_point(aes(y = all, color = "unstructured data")) +  # Plot for data
    labs(title = "", x = "Time", y = "Total Cases") +
    geom_line(data = solns_one, aes(x = time, y = X, color="fitted")) +  # Fitted line for nonstructured cases
    geom_line(data = projected_solns_one_plot, aes(x = time, y = X, color="projected"), linetype="dashed") +  # projection
    scale_y_log10() +
    geom_vline(xintercept=t_split + 0.1, linetype="dotted", color="gray22", linewidth = 0.7) +
    scale_color_manual(name = "Legend",
                       values = c("unstructured data" = "gray22",
                                  "fitted" = "gray22",
                                  "projected" = "gray22")) +
    theme_minimal()
  
  grid.arrange(fig + xlim(0,30) + ylim(0,500), 
               fig_un + xlim(0,30) + ylim(0,500), ncol=1)
  # -------------------------------------------------------------------------- #
  
  MAEadult <- mean(abs(struct_split$testing_data$adults - projected_solns$adults))
  MAEchild <- mean(abs(struct_split$testing_data$children - projected_solns$children))

  mae <- mean(abs((struct_split$testing_data$adults + struct_split$testing_data$children) -
                    (projected_solns$adults + projected_solns$children)))

  mae_un <- mean(abs(unstruct_split$testing_data$all - projected_solns_one$X))

  return(list(fig=fig, fig_un=fig_un, mae=mae, mae_un=mae_un))
}
# ============================================================================ #