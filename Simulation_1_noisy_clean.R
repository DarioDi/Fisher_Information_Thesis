#Simulation 1 - e is changing parameter and time step is 600 (model_parameters_1) ####

library(deSolve)
library(ggplot2)
library(rootSolve)
library(mgcv)
library(dplyr)
library(tibbletime)

source("helper_functions.R")
source("model_parameters_1.R")

#set seed to ensure reproducibility

set.seed(10)

#set parameters for changing observation error and rate of change

sim1_parameter = expand.grid(obs_error = seq(0.00, 0.3, length.out = 7),
                             rate_of_change = seq(0.0006, 0.0016, length.out = 11),
                             replicate = 1:10)

sim1_parameter$replicate <- 1:nrow(sim1_parameter)

n_sim1 = nrow(sim1_parameter)

for(i in seq_len(length.out = n_sim1)){
  
  #run the simulation
  
  sim_current <- run_simulation(parameters = sim1_parameter[i,],
                                time_seq = time_series)
  
  #fit population density of the simulation using GAM
  
  sim_predictions <- extract_gam_predictions(parameters = sim1_parameter[i,], 
                                             sim = sim_current,
                                             time_seq = time_series)
  
  #calculate fisher information of the system
  
  sim_fisher_current <- calc_fisher_current(parameters = sim1_parameter[i,],
                                            predictions = sim_predictions,
                                            time_seq = time_series)
  
  #calculate fisher information of the system using log transform of the s
  #ystem population densities
  
  sim_log_fisher_current <- calc_fisher_current(parameters = sim1_parameter[i,],
                                                predictions = log(sim_predictions),
                                                time_seq = time_series)
  
  #find the min fisher information timepoint for when fisher information d
  #etects the regime shift
  
  sim_fisher_min <- which.min(sim_fisher_current)
  sim1_parameter[i, "sim_fisher_min"] = sim_fisher_min
  
  #find the min fisher information timepoint for when fisher information 
  #of the log transform data detects the regime shift
  
  sim_log_fisher_min <- which.min(sim_log_fisher_current)
  sim1_parameter[i, "sim_log_fisher_min"] = sim_log_fisher_min
  
  #find the regime shift of the simulation
  
  regime_shift_time <- calc_regime_shift(parameters = sim1_parameter[i,],
                                         times = time_series)
  sim1_parameter[i, "regime_shift_time"] = regime_shift_time
  
}

#table for simulation data collected

View(sim1_parameter)

#set up the summary of the data so that the mean for fisher information,
#log fisher information and regime shift are calculated from the replicates
#of the simulations

sim1_parameter_summary = sim1_parameter %>%
  group_by(obs_error, rate_of_change)%>%
  summarize(fisher_min_mean = mean(sim_fisher_min),
            log_fisher_min_mean = mean(sim_log_fisher_min),
            regime_shift_time_mean = mean(regime_shift_time),
            fisher_diff_mean = mean(regime_shift_time-sim_fisher_min),
            fisher_log_diff_mean = mean(regime_shift_time-sim_log_fisher_min))

#plot the difference between fisher information and the regime shift

ggplot(sim1_parameter_summary, aes(y = obs_error, x = rate_of_change, 
                                   fill = fisher_diff_mean)) +
  ggtitle("Fisher Information Difference, Time step of 600")+
  xlab("rate of change of parameter e")+
  ylab("observation error")+
  geom_tile()+
  scale_fill_gradient2()

#plot the difference between the log fisher information and the regime shift

ggplot(sim1_parameter_summary, aes(y = obs_error, x = rate_of_change,
                                   fill = fisher_log_diff_mean)) +
  ggtitle("Log Fisher Information Difference, Time step of 600")+
  xlab("rate of change of parameter e")+
  ylab("observation error")+
  geom_tile()+
  scale_fill_gradient2()




