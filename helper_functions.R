
# Primary workflow functions ----------------------------------------------

#Running the simulation ####

run_simulation <- function(parameters,   # single row dataframe with the varying parameters
                           time_seq, 
                           init_cond = init_cond_default, # Init conditions
                           other_params = params_unchanging){
  
  # 1. Run and the solve equation system
  # Create rate of change timeseries
  rate_1 <- parameters$rate_of_change
  
  # Run the solver
  sim = deSolve::ode(y=init_cond,
                     times = time_seq,
                     func = troph_tri_static_1,
                     parms = append(other_params, 
                                    list(rate_from_time = rate_from_time, 
                                         rate_1 = rate_1)))
  
  # 2. return a data.frame with the 3 columns P, F, J
  return(sim)
}

#Fitting population densities of the trophic triangle using GAM ####

extract_gam_predictions <- function(parameters, sim, time_seq){
  
  
  
  
  # 1. apply the observation error
  err <- parameters$obs_error
  
  n_steps <- length(time_seq)
  
  sim[,"P"] = sim[,"P"]*exp(rnorm(n_steps,0,sd = err))
  sim[,"J"] = sim[,"J"]*exp(rnorm(n_steps,0,sd = err))
  sim[,"F"] = sim[,"F"]*exp(rnorm(n_steps,0,sd = err))
  
  # 2. Run gams
  
  gam_P = gam(P~s(time_seq, k= 20, bs= "ad"),data= as.data.frame(sim),
              family = Gamma(link = "log"),method="REML")
  
  gam_J = gam(J~s(time_seq, k= 20, bs= "ad"),data= as.data.frame(sim),
              family = Gamma(link = "log"),method="REML")
  
  gam_F = gam(F~s(time_seq, k= 20, bs= "ad"),data= as.data.frame(sim),
              family = Gamma(link = "log"),method="REML")
  
  # 3. Extract
  
  predictions <- data.frame(P = fitted(gam_P),
                            F = fitted(gam_F), 
                            J = fitted(gam_J))
  
  return(predictions)
}

#Calculating Fisher Information ####

calc_fisher_current <- function(parameters, predictions, time_seq, 
                                first_deriv_func = calc_1st_deriv, 
                                second_deriv_func = calc_2nd_deriv){
  
  n_steps <- length(time_seq)
  
  #1. Fill first and second derivative matrices
  
  delta <-time_seq[2] - time_seq[1]
  
  first_deriv <- matrix(NA, nrow = n_steps, ncol =3)
  first_deriv[,1] <- first_deriv_func(predictions[,"P"], delta)
  first_deriv[,2] <- first_deriv_func(predictions[,"F"], delta)
  first_deriv[,3] <- first_deriv_func(predictions[,"J"], delta)
  
  second_deriv <- matrix(NA, nrow = n_steps, ncol =3)
  second_deriv[,1] <- second_deriv_func(predictions[,"P"], delta)
  second_deriv[,2] <- second_deriv_func(predictions[,"F"], delta)
  second_deriv[,3] <- second_deriv_func(predictions[,"J"], delta)
  
  fisher_info <- matrix(NA, nrow = n_steps, ncol=1)
  
  #2. set up calculation of fisher information
  
  for(i in seq_len(n_steps)){
    numerator <-  sum(first_deriv[i,]*second_deriv[i,])^2
    denominator <- sqrt(sum(second_deriv[i,]^2))^6
    fisher_info[i,] <-  numerator/denominator
  }
  
  #3. use rolling mean on fisher info and extract
  
  rolling_mean <- rollify(mean, window = 20)
  rolling_mean_fisher <- rolling_mean(fisher_info[,1])
  
  return(rolling_mean_fisher)
  
}


#Determining when the Regime Shift Occurs ####




calc_regime_shift <- function(parameters, times,
                              other_params = params_unchanging,
                              init_cond = init_cond_default){
  
  rate_1 <- parameters$rate_of_change
  model_parameters_1 <- append(other_params, 
                               list(rate_from_time = rate_from_time, 
                                    rate_1 = rate_1))
  
  
  #create a vector of 600 time steps
  n_steps <- length(times)
  
  #create a data frame to hold equilibrium values
  stable_states_1 <- data.frame(P = rep(0,times = n_steps),
                                F = rep(0,times = n_steps),
                                J = rep(0,times = n_steps),
                                eigen  = rep(0,times = n_steps),
                                time = times)
  
  #set the current state for the loop to the original initial condition
  current_state_1 <- init_cond
  
  for(i in seq_len(n_steps)){
    current_time <- times[i]
    #calculate the closest equilibrium point at the current time step
    root_value_1 <- stode(y= current_state_1,
                          time =current_time,
                          func = troph_tri_static_1,
                          jacfunc = troph_tri_jacobian_1,
                          parms = model_parameters_1,
                          positive = TRUE #this ensures that rootSolve will only find positive (or zero) solutions
    )
    
    #change the current state to this value
    current_state_1 <- root_value_1$y
    stable_states_1[i, c("P","F","J")] <- current_state_1
    
    #Calculate the Jacobian of the system at this equilibrium
    current_jacobian_1 <- troph_tri_jacobian_1(t = current_time, 
                                               y = current_state_1,
                                               parms = model_parameters_1)
    
    #calculate eigenvalues of this Jacobian and find the maximum real eigenvalue
    current_eigs_1 <- eigen(current_jacobian_1)
    stable_states_1[i,"eigen"] <- max(Re(current_eigs_1$values))
    
    #add a small perturbation to the current state to keep rootSolve from finding
    #only zero values after the regime shift.
    current_state_1 <- current_state_1 +0.5
  }
  
  #Find the regime shift point as the place where the eigen value of the Jacobian goes to zero (or just above)
  regime_shift <- stable_states_1$time[stable_states_1$eigen==max(stable_states_1$eigen)]
  if(length(regime_shift)>1) warning("multiple regime shifts found") 
  
  return(min(regime_shift)) # TODO verify if minimum is appropriate here
  
}


# Secondary ---------------------------------------------------------------


#Calculating Fisher Information ####

calc_1st_deriv = function(y,delta){
  (lead(y,1) - lag(y,1))/(2*delta)
}

calc_2nd_deriv = function(y,delta){
  (lead(y,1) + lag(y,1)-2*y)/delta^2
}

