
# Primary workflow functions ----------------------------------------------

#Running the simulation ####

params_unchanging <- 
  list(T_mat = 5,    #length of time it takes a juvenile predator to mature to an adult
       m = 0.025,    #mortality rate of adult and juvenile fish
       s = 0.05,     #The stocking rate for adult predators
       
       f = 0.5,      #amount of new offspring for each adult predator per unit time
       a_PJ  = 0.05, #Cannibalism rate of adult predators on juveniles
       a_FJ  = 0.1,  #attack rate of forage fish on juvenile predators
       
       r = 0.25,     #population growth rate of forage fish at low densities
       b = 0.005,    #density-dependence term for the forage fish
       a_PF = 0.1,   #attack rate of adult predators on forage fish
       d = 0.5       #Stocking rate for forage fish
  ) 

time_series <- seq(1,100,length.out = 100)
init_cond_default <- c(P = 77, F = 0.067, J = 9.37)


run_simulation <- function(parameters,   # single row dataframe with the varying parameters
                           time_seq = time_series, 
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

extract_gam_predictions <- function(parameters, sim, time_seq = time_series){
  
  set.seed(10)
  
  
  # 1. apply the observation error
  err <- parameters$obs_error
  
  sim[,"P"] = sim[,"P"]*exp(rnorm(100,0,sd = err))
  sim[,"J"] = sim[,"J"]*exp(rnorm(100,0,sd = err))
  sim[,"F"] = sim[,"F"]*exp(rnorm(100,0,sd = err))
  
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

calc_fisher_current <- function(parameters, predictions, 
                                time_seq = time_series, 
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
  
  rolling_mean <- rollify(mean, window = 10)
  rolling_mean_fisher <- rolling_mean(fisher_info[,1])
  
  return(rolling_mean_fisher)
  
}


#Determining when the Regime Shift Occurs ####




calc_regime_shift <- function(parameters, other_params = params_unchanging,
                              init_cond = init_cond_default, times  = time_series){
  
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
  
  return(min(regime_shift)) # TODO verify if minimum is appropriate here
  
}


# Secondary ---------------------------------------------------------------

#Running the simulation ####

troph_tri_static_1 = function(t,y,parms){
  
  #This code extracts the three state variables (P,F, and J) from the y vector
  #so they can be referred to as single letters in the equations for rates of
  #change
  P = y["P"]
  F = y["F"]
  J = y["J"]
  
  #This next code calculates the derivatives at each point in time. 
  #the with(x,...) function here make the model parameters available by name
  #without having to type parms$e*P + parms$s...
  dP = with(parms, J/T_mat - m*P - rate_from_time(t, rate_1)*P + s) 
  dF = with(parms, r*F - b*F^2 - a_PF*P*F + d)
  dJ = with(parms, f*P - J/T_mat - m*J - a_PJ*P*J - a_FJ*F*J)
  return(list(c(P=dP,F=dF, J=dJ)))
}

rate_from_time <- function(rate_1, t, 
                           min_value_1=0, max_value_1=1, lag_time_1=0){
  value_1=rate_1*t-lag_time_1
  value_1[value_1 < min_value_1] <- min_value_1
  value_1[value_1 > max_value_1] <- max_value_1
  value_1
}

#Calculating Fisher Information ####

calc_1st_deriv = function(y,delta){
  (lead(y,1) - lag(y,1))/(2*delta)
}

calc_2nd_deriv = function(y,delta){
  (lead(y,1) + lag(y,1)-2*y)/delta^2
}

#Calculating Regime Shift point ####

troph_tri_jacobian_1 <- function(t,y,parms){
  
  jacobian_1 <- matrix(NA, nrow=3, ncol=3)
  
  #dP - differentiation variable: P (-e(t)-m) 
  jacobian_1[1,1] <- with(parms, (-rate_from_time(t, rate_1)-m))
  #dP - differentiation variable: F 0
  jacobian_1[1,2] <- (0)
  #dP - differentiation variable: J (1/T_mat)
  jacobian_1[1,3] <- with(parms, (1/T_mat))
  
  
  #dF - differentiation variable: P (-a_PF*F)
  jacobian_1[2,1] <- with(parms, (-a_PF*y[2]))
  #dF - differentiation variable: F (-2*b*F)+(r)-(a_PF*P)
  jacobian_1[2,2] <- with(parms, ((-2*b*y[2])+(r)-(a_PF*y[1])))
  #dF - differentiation variable: J 0
  jacobian_1[2,3] <- (0)
  
  #dJ - differentiation variable: P (F-(a_PJ*J)
  jacobian_1[3,1] <- with(parms, (f-(a_PJ*y[3])))
  #dJ - differentiation variable: F (P-(a_FJ*J)
  jacobian_1[3,2] <- with(parms, -a_FJ*y[3])
  #dJ - differentiation variable: J (-1/T_mat)-(a_PJ*P)-(m)-(a_FJ*F)
  jacobian_1[3,3] <- with(parms, ((-1/T_mat)-(a_PJ*y[1])-(m)-(a_FJ*y[2])))
  
  return(jacobian_1)
}

