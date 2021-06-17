# Primary workflow functions ----------------------------------------------

#set time step of the simulation ####

time_series <- seq(1,600,length.out = 600)

#Running the simulation ####

params_unchanging <- 
  list(T_mat = 5,    #length of time it takes a juvenile predator to mature to an adult
       m = 0.025,    #mortality rate of adult and juvenile fish
       s = 0.05,     #The stocking rate for adult predators
       e = 0.005,    #The extraction rate at any given time
       
       f = 0.5,      #amount of new offspring for each adult predator per unit time
       a_PJ  = 0.05, #Cannibalism rate of adult predators on juveniles
       a_FJ  = 0.1,  #attack rate of forage fish on juvenile predators
       
       r = 0.25,     #population growth rate of forage fish at low densities
       b = 0.005,    #density-dependence term for the forage fish
       d = 0.5       #Stocking rate for forage fish
  ) 

init_cond_default <- c(P = 77, F = 0.067, J = 9.37)

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
  dP = with(parms, J/T_mat - m*P - e*P + s) 
  dF = with(parms, r*F - b*F^2 - rate_from_time(t, rate_1)*P*F + d)
  dJ = with(parms, f*P - J/T_mat - m*J - a_PJ*P*J - a_FJ*F*J)
  return(list(c(P=dP,F=dF, J=dJ)))
}

rate_from_time <- function(rate_1, t, 
                           min_value_1=0, max_value_1=1, lag_time_1=0.1){
  value_1=lag_time_1+rate_1*t
  value_1[value_1 < min_value_1] <- min_value_1
  value_1[value_1 > max_value_1] <- max_value_1
  value_1
}

#Calculating Regime Shift point ####

troph_tri_jacobian_1 <- function(t,y,parms){
  
  jacobian_1 <- matrix(NA, nrow=3, ncol=3)
  
  #dP - differentiation variable: P 
  jacobian_1[1,1] <- with(parms, (-m-e))
  #dP - differentiation variable: F 
  jacobian_1[1,2] <- (0)
  #dP - differentiation variable: J 
  jacobian_1[1,3] <- with(parms, (1/T_mat))
  
  
  #dF - differentiation variable: P 
  jacobian_1[2,1] <- with(parms, (-rate_from_time(t, rate_1)*y[2]))
  #dF - differentiation variable: F - error
  jacobian_1[2,2] <- with(parms, ((-2*b*y[2])-rate_from_time(t, rate_1)*y[1])+(r))
  #dF - differentiation variable: J
  jacobian_1[2,3] <- (0)
  
  
  #dJ - differentiation variable: P 
  jacobian_1[3,1] <- with(parms, (f-(a_PJ*y[3])))
  #dJ - differentiation variable: F 
  jacobian_1[3,2] <- with(parms, (-(a_FJ*y[3])))
  #dJ - differentiation variable: J 
  jacobian_1[3,3] <- with(parms, ((-1/T_mat)-(a_PJ*y[1])-(m)-(a_FJ*y[2])))
  
  return(jacobian_1)
}




