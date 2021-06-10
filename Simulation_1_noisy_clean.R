library(deSolve)
library(ggplot2)
library(rootSolve)
library(mgcv)
library(dplyr)
library(tibbletime)

source("helper_functions.R")

sim1_parameter = expand.grid(obs_error = seq(0.01, 1, length.out = 5),
                             rate_of_change = seq(0.0001, 0.001, length.out = 5))
sim1_parameter$replicate <- 1:nrow(sim1_parameter)

n_sim1 = nrow(sim1_parameter)

for(i in seq_len(length.out = n_sim1)){
  
  # feeding single row df to function
  sim_current <- run_simulation(parameters = sim1_parameter[i,])
  
  print(head(sim_current))
  
  sim_predictions <- extract_gam_predictions(parameters = sim1_parameter[i,], 
                                             sim = sim_current)
  
  print(head(sim_predictions))  
  
  sim_fisher_current <- calc_fisher_current(parameters = sim1_parameter[i,],
                                            predictions = sim_predictions)
  
  # sim1_fisherinfo_current = rolling_mean_fisher_1 #adjust function
  # sim1_fisherinfo_log_current = #create function for FI with log
  #   sim1_fisherinfo_min = 
  #   sim1_fisherinfo_log_min =
  #   sim1_regimeshift_time = regime_shift_1 #need to fix function
  # sim1_parameter[i, "sim1_fisherinfo_min"] = sim1_fisherinfo_min
  # sim1_parameter[i, "sim1_fisherinfo_log_min"] = sim1_fisherinfo_log_min
  # sim1_parameter[i, "sim1_regimeshift_time"] = sim1_regimeshift_time
}






