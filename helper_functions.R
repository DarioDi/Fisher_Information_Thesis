
# Primary workflow functions ----------------------------------------------

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

time_series <- seq(0,600,length.out = 100)


run_simulation <- function(parameters,   # single row dataframe with the varying parameters
                           time_seq = time_series, 
                           init_cond = c(P = 77, F = 0.067, J = 9.37), # Init conditions
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

extract_gam_predictions <- function(parameters, sim, time_seq = time_series){
  
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

# Secondary ---------------------------------------------------------------

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