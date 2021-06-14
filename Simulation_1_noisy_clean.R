library(deSolve)
library(ggplot2)
library(rootSolve)
library(mgcv)
library(dplyr)
library(tibbletime)

source("helper_functions.R")

sim1_parameter = expand.grid(obs_error = seq(0.00, 0.3, length.out = 7),
                             rate_of_change = seq(0.0005, 0.0015, length.out = 11))

sim1_parameter$replicate <- 1:nrow(sim1_parameter)

n_sim1 = nrow(sim1_parameter)

for(i in seq_len(length.out = n_sim1)){
  
  # feeding single row df to function
  sim_current <- run_simulation(parameters = sim1_parameter[i,])
  
 # print(head(sim_current))
  
  sim_predictions <- extract_gam_predictions(parameters = sim1_parameter[i,], 
                                             sim = sim_current)
  
 # print(head(sim_predictions))  
  
  sim_fisher_current <- calc_fisher_current(parameters = sim1_parameter[i,],
                                            predictions = sim_predictions)
  
  #print(head(sim_fisher_current))
  
  sim_log_fisher_current <- calc_fisher_current(parameters = sim1_parameter[i,],
                                                predictions = log(sim_predictions))
  
  #print(head(sim_log_fisher_current))
  
  sim_fisher_min <- which.min(sim_fisher_current)
  sim1_parameter[i, "sim_fisher_min"] = sim_fisher_min
  
  sim_log_fisher_min <- which.min(sim_log_fisher_current)
  sim1_parameter[i, "sim_log_fisher_min"] = sim_log_fisher_min
  
  regime_shift_time <- calc_regime_shift(parameters = sim1_parameter[i,])
  sim1_parameter[i, "regime_shift_time"] = regime_shift_time
  
}

View(sim1_parameter)


ggplot(sim1_parameter, aes(y = obs_error, x = rate_of_change, fill = sim_fisher_min)) +
  geom_tile()
ggplot(sim1_parameter, aes(y = obs_error, x = rate_of_change, fill = sim_log_fisher_min)) +
  geom_tile()
ggplot(sim1_parameter, aes(y = obs_error, x = rate_of_change, fill = regime_shift_time)) +
  geom_tile()

