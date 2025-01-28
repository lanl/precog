# generate synthetic phylogenies and time series
# call script plot_single_projection.R (splits data into training/testing,
# plots the fits and projections, calculates MAEs)
# fit time series models to training data and estimate params from time series

# code by Julie A. Spencer and Scottie Alexander
# 12/11/2024

library(TreeSim)
library(ape)
library(geiger)
library(devtools)
library(remotes)
library(TreePar)
library(phytools)
library(bbmle)
library(stats4)
library(fitode)
library(deSolve)
library(Deriv)
library(plyr)
library(tidyverse)
library(gridExtra)

source("./plot_single_projection.R")

# simulate additional trees for time series fitting and projection
# estimate parameters on training data and test on holdout data
# first with age-structured ODE model, then with homogeneous ODE model

# global params
set.seed(42)

n<-1000
numbsim <- 200
init<- -1 
sampprob1<-1
sampprob2<-1
s<-c(sampprob1,sampprob2)

### Sc1 param set-- homogeneous transmission ###
lambda11_H0<-3
lambda12_H0<-3
lambda21_H0<-3
lambda22_H0<-3
death1_H0<-3
death2_H0<-3
l_H0<-rbind(c(lambda11_H0,lambda12_H0),c(lambda21_H0,lambda22_H0))
d_H0<-c(death1_H0,death2_H0)

#Sc2 param set--transmission within groups >> transmission between groups ###
lambda11_H1<-5
lambda12_H1<-1
lambda21_H1<-1
lambda22_H1<-5
death1_H1<-3
death2_H1<-3
l_H1<-rbind(c(lambda11_H1,lambda12_H1),c(lambda21_H1,lambda22_H1))
d_H1<-c(death1_H1,death2_H1)
     
### Sc3 param set -- asymmetric: transmission high from adults, low from children###
lambda11_H2<-5
lambda12_H2<-5
lambda21_H2<-1
lambda22_H2<-1
death1_H2<-3
death2_H2<-3
l_H2<-rbind(c(lambda11_H2,lambda12_H2),c(lambda21_H2,lambda22_H2))
d_H2<-c(death1_H2,death2_H2)

### Sc4 param set -- extra asymmetric: transmission high only within adults;
### low within children and among groups ###
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

# simulate 200 trees with 1000 tips for each scenario
trees_H0pr<-lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
                 lambdavector=l_H0, deathvector=d_H0, sampprobvector=s,
                 init=init)
trees_H1pr<-lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
                 lambdavector=l_H1, deathvector=d_H1, sampprobvector=s,
                 init=init)
trees_H2pr<-lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
                 lambdavector=l_H2, deathvector=d_H2, sampprobvector=s,
                 init=init)
trees_H3pr <- lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
                   lambdavector=l_H3, deathvector=d_H3, sampprobvector=s,
                   init=init)

# add root.dist (sum of branch lengths leading to each tip for multiple trees)
trees_H0pr <- lapply(trees_H0pr, add_root_distance)
trees_H1pr <- lapply(trees_H1pr, add_root_distance)
trees_H2pr <- lapply(trees_H2pr, add_root_distance)
trees_H3pr <- lapply(trees_H3pr, add_root_distance)

# function to derive time series from trees
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

# generate 200 time series for each scenario
df_H0pr <- lapply(trees_H0pr, gen_ts)
df_H1pr <- lapply(trees_H1pr, gen_ts)
df_H2pr <- lapply(trees_H2pr, gen_ts)
df_H3pr <- lapply(trees_H3pr, gen_ts)



# for the homogeneous (unstructured) time series, combine age groups and delete age columns
combine <- function(timeseries){
  df <- timeseries %>% mutate(all=adults + children) 
  df <- df[-c(2,3)]
  return(df)
}

df_H0_one <- lapply(df_H0pr, combine)
df_H1_one <- lapply(df_H1pr, combine)
df_H2_one <- lapply(df_H2pr, combine)
df_H3_one <- lapply(df_H3pr, combine)


# for the structured data, combine age groups in a new column, "all."

combine_pr <- function(timeseries){
  df <- timeseries %>% mutate(all=adults + children)
  return(df)
}

df_H0pr <- lapply(df_H0pr, combine_pr)
df_H1pr <- lapply(df_H1pr, combine_pr)
df_H2pr <- lapply(df_H2pr, combine_pr)
df_H3pr <- lapply(df_H3pr, combine_pr)

# ============================================================================ #
### define exponential growth two-state ODE model ###
exp_model_2 <- odemodel(
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
# ============================================================================ #
### define exponential growth one-state ODE model ###
exp_model_1 <- odemodel(
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

#### Illustrative Case Projections Start Here #################################
###############################################################################
#### Projections ###############################################################

# H0 (Scenario 1) TWO-STATE and ONE-STATE -- estimate params and simulate projections
start0 <- c(beta11=.3, beta12=.3, beta21=.3, beta22=.3, I10=1, I20=1,
           size1=3, size2=3, gamma=.3)

start_un0 <- c(r = 0.2, X0=1)

N <- length(df_H0pr)
df_out <- data.frame(
  scenario1_mae = rep(NaN, N),
  scenario1_mae_un = rep(NaN, N),
  scenario2_mae = rep(NaN, N),
  scenario2_mae_un = rep(NaN, N),
  scenario3_mae = rep(NaN, N),
  scenario3_mae_un = rep(NaN, N),
  scenario4_mae = rep(NaN, N),
  scenario4_mae_un = rep(NaN, N)
)

pdf(file = "projection_examples_scenario_1.pdf", width = 6, height = 9)

for (k in 1:length(df_H0pr)) {
  h <- plot_single_projection(df_H0pr[[k]], start0, exp_model_2,
                              df_H0_one[[k]], start_un0, exp_model_1)
  df_out[k,"scenario1_mae"] = h$mae
  df_out[k,"scenario1_mae_un"] = h$mae_un
}
dev.off()

# H1 (Scenario 2) TWO-STATE and ONE-STATE -- estimate params and simulate projections
start1 <- c(beta11=.5, beta12=.1, beta21=.1, beta22=.5, I10=1, I20=1,
           size1=3, size2=3, gamma=.3)

start_un1 <- c(r = 0.2, X0=1)

pdf(file = "projection_examples_scenario_2.pdf", width = 6, height = 9)

for (k in 1:length(df_H1pr)) {
  h <- plot_single_projection(df_H1pr[[k]], start1, exp_model_2,
                              df_H1_one[[k]], start_un1, exp_model_1)
  df_out[k,"scenario2_mae"] = h$mae
  df_out[k,"scenario2_mae_un"] = h$mae_un
}
dev.off()


# H2 (Scenario 3) TWO-STATE and ONE-STATE -- estimate params and simulate projections
start2 <- c(beta11=.5, beta12=.5, beta21=.1, beta22=.1, I10=1, I20=1,
           size1=3, size2=3, gamma=.3)

start_un2 <- c(r = 0.2, X0=1)

pdf(file = "projection_examples_scenario_3.pdf", width = 6, height = 9)

for (k in 1:length(df_H2pr)) {
  h <- plot_single_projection(df_H2pr[[k]], start2, exp_model_2,
                              df_H2_one[[k]], start_un2, exp_model_1)
  df_out[k,"scenario3_mae"] = h$mae
  df_out[k,"scenario3_mae_un"] = h$mae_un
}
dev.off()


# H3 (Scenario 4) TWO-STATE and ONE-STATE--estimate params and project simulations
start3 <- c(beta11=.5, beta12=.1, beta21=.1, beta22=.1, I10=1, I20=1,
             size1=3, size2=3, gamma=.3)

start_un3 <- c(r = 0.2, X0=1)

pdf(file = "projection_examples_scenario_4.pdf", width = 6, height = 9)

for (k in 1:length(df_H3pr)) {
  h <- plot_single_projection(df_H3pr[[k]], start3, exp_model_2,
                              df_H3_one[[k]], start_un3, exp_model_1)
  df_out[k,"scenario4_mae"] = h$mae
  df_out[k,"scenario4_mae_un"] = h$mae_un
}
dev.off()

################################################################################
write.csv(df_out, "results_MAE/results_MAE_12Jan2025")

# plotting code for MAE boxplots

# plot for new MAEs December 9, 2024 (code debugged)
MAE <- read.csv("results_MAE/results_MAE_12Jan2025")

dfH0 <- data.frame(
  "Scenario 1" = MAE$scenario1_mae,
  "Scenario 2" = MAE$scenario2_mae,
  "Scenario 3" = MAE$scenario3_mae,
  "Scenario 4" = MAE$scenario4_mae
)

dfH0un <- data.frame(
  "Scenario 1" = MAE$scenario1_mae_un,
  "Scenario 2" = MAE$scenario2_mae_un,
  "Scenario 3" = MAE$scenario3_mae_un,
  "Scenario 4" = MAE$scenario4_mae_un
)

# format the data for plotting
dfH0_long <- dfH0 %>% pivot_longer(cols = everything(), names_to = "Scenario",
                                   values_to = "MAE")
df_structured <- dfH0_long %>% mutate(Model = "structured")
dfH0un_long <- dfH0un %>% pivot_longer(cols = everything(), names_to = "Scenario",
                                       values_to = "MAE")
df_unstructured <- dfH0un_long %>% mutate(Model = "homogeneous")
all_MAE <- rbind(df_structured,df_unstructured)


# boxplots for MAE distribution, both models
pMAE <- ggplot(all_MAE, aes(x= Scenario, y = MAE, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             size = .4, alpha = 0.4) +  
  labs(title = "Mean Absolute Error of Case Projections") +
  ylim(0,150) +
  scale_fill_brewer(palette= "Paired") +
  theme_minimal()
pMAE

