# minimum reproducible example
# December 6, 2024



# this code will 
# 1. simulate one data set for Scenario 4 (extra asymmetrical),
#    consisting of three parts (each with 1000 tips/cases)
#    (a) an age-structured tree, 
#    (b) an age structured time series
# 2. define the ODE models (structured and nonstructured)
# 3. estimate 4 transmission parameters from the structured trees
#    and time series
# 4. demonstrate the effect of age structure on case projections
#       (i)  create unstructured time series and split 
#               both structured and unstructured sets into training/testing sets
#       (ii) fit models to time series training data (structured and unstructured)
#       (iii) project simulations through the testing data (structured & unstructured)
#       (iv) plot the fits and projections (structured & unstructured)
#       (v) calculate the MAEs for the projections (structured & unstructured)

### load libraries ###
library(TreeSim)
library(ape)
library(geiger)
library(devtools)
library(remotes)
library(TreePar)
library(phytools)
library(tidyverse)
library(bbmle)
library(stats4)
library(fitode)
library(deSolve)
library(Deriv)
library(plyr)
library(gridExtra)

#==============================================================================#
### SECTION ONE--simulate the data: one tree and one time series, structured ###
#==============================================================================#
### define global parameters for synthetic trees and time series ###
set.seed(42)
n<-1000
numbsim <- 1
init<- -1 
sampprob1<-1
sampprob2<-1
s<-c(sampprob1,sampprob2)

### Scenario 4 param set--extra asymmetric ###
# transmission within ONE GROUP >> 3 other rates #
# extra asymmetric #
lambda11_H3<-5
lambda12_H3<-1
lambda21_H3<-1
lambda22_H3<-1
death1_H3<-3
death2_H3<-3
l_H3<-rbind(c(lambda11_H3,lambda12_H3),c(lambda21_H3,lambda22_H3))
d_H3<-c(death1_H3,death2_H3)

# root distance function (for time of sample)
add_root_distance <- function(phy)
{
  # order = nodes out from root, then tips
  i_node <- c((Ntip(phy)+1) : (Ntip(phy)+Nnode(phy)), seq(Ntip(phy)))
  phy$root.dist <- rep(NA, length(i_node))
  
  # at the root
  phy$root.dist[i_node[1]] <- 0
  
  # all other nodes and tips
  for (n in i_node[-1]) {
    phy$root.dist[n] <- phy$edge.length[which(phy$edge[,2] == n)] + 
      phy$root.dist[getParent(phy, n)]
  }
  
  return(phy)
}

# simulate one tree
tree_one <- lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
                   lambdavector=l_H3, deathvector=d_H3, sampprobvector=s,
                   init=init)

# add root distance to one tree
tree_one <- lapply(tree_one, add_root_distance)
# note--the phylo object will be in the global environment, not saved to disk.

# function to count cases (tips from the trees generated above)
# and time-bin the cases
gen_ts <- function(tree){
  ts <- data.frame(tree$tip.label,
                   tree$root.dist[1:Ntip(tree)],
                   tree$state)
  colnames(ts) <- c("tip", "time", "state")
  ts <- ts %>% mutate(cases = 1)
  ts <- pivot_wider(ts, names_from ="state", values_from ="cases")
  ts[is.na(ts)] <- 0
  colnames(ts) <- c("tip", "time", "adults", "children")
  ts$time <- ts$time*10
  ts$time <- as.integer(ts$time)
  ts$adults <- as.integer(ts$adults)
  ts$children <- as.integer(ts$children)
  ts <- ts[-1]
  ts <- plyr::ddply(ts,"time",numcolwise(sum))
  return(ts)
}

# call function to generate one time series
df_one <- lapply(tree_one, gen_ts) # creates a list 

# check for 1000 tips
num_tips <- sum(df_one[[1]][["adults"]]) + 
  sum(df_one[[1]][["children"]])
num_tips
# check data types
glimpse(df_one)

# save stuctured time series as csv file
write.csv(df_one,
            paste0('ts_files/df_one.csv'),
            row.names=FALSE)

#### end of data simulation section ####

#===========================================================================#
### SECTION TWO_1: define the ODE models ###
#===========================================================================#
### define exponential growth two-state ODE model ###
exp_2state <- odemodel(
  name="Two-State Transmission Model",
  model=list(
    dI1 ~ I1*beta11 + I2*beta21 - I1*gamma,
    dI2 ~ I2*beta22 + I1*beta12 - I2*gamma
  ),
  observation=list(
    adults ~ dnbinom(mu=I1, size=size1),
    children ~ dnbinom(mu=I2, size=size2)
  ),
  initial=list(
    I1 ~ I10,
    I2 ~ I20
  ),
  par=c("beta11", "beta12", "beta21", "beta22",
        "I10", "I20", "size1", "size2", "gamma")
)

exp_1state <- odemodel(
  name="exponential",
  model=list(
    X ~ r * X
  ),
  observation=list(
    all ~ dpois(lambda=X)
  ),
  initial=list(
    X ~ X0
  ),
  par=c("r", "X0")
)

### SECTION TWO_2: define functions that evaluate the ODE models ###
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

### define exponential growth one-state ODE model ###
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
#### end of model definition and function section ####

#===========================================================================#
### SECTION THREE: estimate parameters from structured trees and time series###
#===========================================================================#

# estimate parameters from one tree (age-structured)
param_one <- matrix(0.0, 1, 5)
for (i in 1){
  p_one <- optim(c(3,3,3,3,3),LikTypesSTT, 
                 phylo=tree_one[[i]],fix=rbind(c(6,7,8),c(-5,0,0),c(1,1,1)),
                 sampfrac=s,survival=0,posR=0,control=list(maxit=10000))
  
  param_one[i,] <- p_one$par
}

result_one <- data.frame(param_one)
write.csv(result_one, "result_one/tree_one.csv", row.names=FALSE)

# estimate parameters from one time series (age-structured)
ts_one <- read.csv("ts_files/df_one.csv")
start_one <- c(beta11=.3, beta21=.3, beta12=.3, beta22=.3, I10=1, I20=1,
            size1=3, size2=3, gamma=.3)

tsfit_one <- fitode(exp_2state, data = ts_one,
                       start = start_one,
                       tcol = "time"
)

ans <- tsfit_one@coef*10
ans
write.csv(ans, "result_one/ts_one.csv")

#### end of parameter estimation section ####

#===========================================================================#
### SECTION FOUR: demonstrate the effect of age structure on case projections
#===========================================================================#

# (i) create unstructured time series 
# and split time series data into training/testing #
# where training set = 0 to 300 cases and testing set = 301-1000 cases.

# to create unstructured time series, combine age groups and delete age columns
# there will now be two columns, "time" and "all" 
combine <- function(timeseries){
  df <- timeseries %>% mutate(all=adults + children) 
  df <- df[-c(2,3)]
  return(df)
}
df_one_unstruct <- lapply(df_one, combine)
# convert to dataframe
df_one_unstruct <- df_one_unstruct[[1]]

# also create "all" column for structured data for finding thresshold
combine_pr <- function(timeseries){
  df <- timeseries %>% mutate(all=adults + children)
  return(df)
}
df_one_struct <- lapply(df_one, combine_pr)
# convert to dataframe
df_one_struct <- df_one_struct[[1]]

# sum cases, find the threshold, and split the data into training and testing

### unstructured train/test split ###
# Create the cumulative sum column 'sum'
df_one_unstruct$sum <- cumsum(df_one_unstruct$all)
# Find the first time step where the cumulative sum is >= threshold
threshold_time <- df_one_unstruct$time[min(which(df_one_unstruct$sum >=300))]
# Create the training dataset (from the first time step to the threshold time step)
training_unstruct <- df_one_unstruct[df_one_unstruct$time < threshold_time, ]
# remove "sum" column to prep for fitting DEBUGGING STEP
training_unstruct <- select(training_unstruct, -sum)
# Create the testing dataset (after the threshold time step)
testing_unstruct <- df_one_unstruct[df_one_unstruct$time >= threshold_time, ]
# remove "sum" column to prep for fitting DEBUGGING STEP
testing_unstruct <- select(testing_unstruct, -sum)

### structured train/test split ###
# Create the cumulative sum column 'sum'
df_one_struct$sum <- cumsum(df_one_struct$all)
# Find the first time step where the cumulative sum is >= threshold
threshold_time <- df_one_struct$time[min(which(df_one_struct$sum >=300))]
# Create the training dataset (from the first time step to the threshold time step)
training_struct <- df_one_struct[df_one_struct$time < threshold_time, ]
# remove "sum" column to prep for fitting DEBUGGING STEP
training_struct <- select(training_struct, -sum)
# remove "all" column to prep for fitting DEBUGGING STEP
training_struct <- select(training_struct, -all)
# Create the testing dataset (after the threshold time step)
testing_struct <- df_one_struct[df_one_struct$time >= threshold_time, ]
# remove "sum" column to prep for fitting DEBUGGING STEP
testing_struct <- select(testing_struct, -sum)
# remove "all" column to prep for fitting DEBUGGING STEP
testing_struct <- select(testing_struct, -all)

# (ii) fit models to time series training data and estimate parameters

##### unstructured #####
# check data format
#view(SierraLeone2014)

# define vector of parameters
start_unstruct <- c(r = 0.2, X0=1)
# fit and estimate 
tsfit_unstruct <- fitode(
  exp_1state, 
  data = subset(training_unstruct,time>=0),
  start = start_unstruct,
  tcol = "time"
)
coef_unstruct <- tsfit_unstruct@coef
coef_unstruct

##### structured #####
# define vector of parameters
start_struct <- c(beta11=.5, beta12=.1, beta21=.1, beta22=.5, I10=1, I20=1,
            size1=3, size2=3, gamma=.3)
tsfit_struct <- fitode(exp_2state, data = training_struct,
                start = start_struct,
                tcol = "time"
)

coef_struct <- tsfit_struct@coef
coef_struct

# (iii) project simulations through the testing data #
### unstructured ###
# Set up a time vector
t1 <- tail(training_unstruct$time, n=1)
time = seq(0, t1, by=1) # in weeks

# Set parameter values
  r = as.numeric(coef_unstruct["r"])
  X0 = as.numeric(coef_unstruct["X0"])
  params_un = c(r, X0)

# Set initial values for state variables
xstart_un = c(X=X0)

# solve
solns_un = as.data.frame(ode(xstart_un, time, exponential_1state, params_un))

# Define the projection time range (for the testing data)
t0 <- testing_unstruct$time[1] - 1
t1 <- tail(testing_unstruct$time, n=1)
projection_times_un <- seq(t0, t1, by = 1)

# Simulate the ODE system as a "future" projection, starting with the last fitted
# value in "solns one"
X0 = tail(solns_un$X, n=1)
xstart_un = c(X=X0)
projected_solns_un <- as.data.frame(ode(xstart_un, projection_times_un,
                                         exponential_1state, params_un))


##### structured #####
# Set up a time vector
t1 <- tail(training_struct$time, n=1)
time_struct = seq(0, t1, by=1) # in weeks

params_struct <- as.numeric(c(coef_struct[1:4], coef_struct["gamma"]))
xstart_struct = c(I1 = as.numeric(coef_struct["I10"]),  I2 = as.numeric(coef_struct["I20"]))
solns = as.data.frame(ode(xstart_struct, time_struct, exponential_2state, params_struct))

# Define the projection time range (the testing data)
t0 <- testing_struct$time[1] - 1
t1 <- tail(testing_struct$time, n=1)
projection_times <- seq(t0, t1, by = 1)

# Simulate the ODE system for projection, starting with the last fitted values
# in "solns"
xstart_p = c(I1 = tail(solns$I1, n=1) , I2 = tail(solns$I2, n=1))
projected_solns <- as.data.frame(ode(xstart_p, projection_times, exponential_2state, params_struct))

# (iv) plot the fits and projections #
### structured ###
fig <- ggplot(df_one_struct, aes(x = time)) +
  geom_point(aes(y = adults, color = "adults-cases")) +  # Plot for data1
  geom_point(aes(y = children, color = "children-cases")) +  # Plot for data2
  labs(title = "", x = "Time", y = "Cases") +
  geom_line(data = solns, aes(x = time, y = I1, color = "adults-fitted")) +  # Fitted line for Data1
  geom_line(data = solns, aes(x = time, y = I2, color = "children-fitted")) +  # Fitted line for Data2
  geom_line(data = projected_solns, aes(x = time, y = I1, color = "adults-projected"),linetype="dashed") + # I1 forecast
  geom_line(data = projected_solns, aes(x = time, y = I2, color = "children-projected"),linetype = "dashed") + #I2 forecast
  geom_vline(xintercept=t0 + 0.1, linetype="dotted", color="gray22", linewidth = 0.7) +
  scale_color_manual(name = "Legend",
                     values = c("adults-cases" = "blue",
                                "children-cases" = "orange",
                                "adults-fitted" = "blue2",
                                "children-fitted" = "darkorange2",
                                "adults-projected" = "blue2",
                                "children-projected" = "darkorange2")) +
  theme_minimal()

fig

### unstructured ###
fig_un <- ggplot(df_one_unstruct, aes(x = time)) +
  geom_point(aes(y = all, color = "unstructured data")) +  # Plot for data
  labs(title = "", x = "Time", y = "Total Cases") +
  geom_line(data = solns_un, aes(x = time, y = X, color="fitted")) +  # Fitted line for nonstructured cases
  geom_line(data = projected_solns_un, aes(x = time, y = X, color="projected"), linetype="dashed") +  # projection
  geom_vline(xintercept=t0 + 0.1, linetype="dotted", color="gray22", linewidth = 0.7) +
  scale_color_manual(name = "Legend",
                     values = c("unstructured data" = "gray22",
                                "fitted" = "gray22",
                                "projected" = "gray22")) +
  theme_minimal()

fig_un

### combine plots ###
grid.arrange(fig + xlim(0,30) + ylim(0,200), 
             fig_un + xlim(0,30) + ylim(0,200), ncol=1)


# (v) calculate the Mean Absolute Error (MAE) #

# for the structured data, combine age groups in a new column, "all." 
# this will be used in the MAE calculation.
# the structure ts will now have 4 columns, "time", "adults", "children", "all".
combine_pr <- function(timeseries){
  df <- timeseries %>% mutate(all=adults + children)
  return(df)
}
df_one_struct <- lapply(df_one, combine_pr)
# convert to dataframe
df_one_struct <- df_one_struct[[1]]

### unstructured ###
# remove first data point of projected solutions (included for plotting purposes)
# so lengths will match and only the actual projected solutions are used in the calculation
projected_solns_un_rem <- projected_solns_un[-1, ]
MAE_un <- mean(abs(testing_unstruct$all - projected_solns_un_rem$X))
MAE_un

### structured ###
# remove first data point of projected solutions (included for plotting purposes)
# so lengths will match and only the actual projected solutions are used in the calculation
projected_solns_rem <- projected_solns[-1, ]
MAEadult <- mean(abs(testing_struct$adults - projected_solns_rem$I1))
MAEchild <- mean(abs(testing_struct$children - projected_solns_rem$I2))
MAEstructured <-  mean(abs((testing_struct$adults + testing_struct$children) -
                                     (projected_solns_rem$I1 + projected_solns_rem$I2)))


