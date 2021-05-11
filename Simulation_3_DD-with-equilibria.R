#Simulation 2 - decreasing the attack rate of predators on forage fish slowly over time


library(deSolve)
library(ggplot2)
library(rootSolve)

#set model parameters

model_parameters_3 = list(
  T_mat = 5,    #length of time it takes a juvenile predator to mature to an adult
  m = 0.025,    #mortality rate of adult and juvenile fish
  s = 0.05,     #The stocking rate for adult predators
  e = function(t,min_value_3=0,max_value_3=1,rate_3=0.001,lag_time_3=0.15){value_3=rate_3*t-lag_time_3
  if(value_3<min_value_3)value_3=min_value_3
  if(value_3>max_value_3)value_3=max_value_3
  value_3
  },      #The extraction rate at any given time
  
  f = 0.5,      #amount of new offspring for each adult predator per unit time
  a_PJ  = 0.05, #Cannibalism rate of adult predators on juveniles
  a_FJ  = 0.1,  #attack rate of forage fish on juvenile predators
  
  r = 0.25,     #population growth rate of forage fish at low densities
  b = 0.005,    #density-dependence term for the forage fish
  a_PF = function(t,min_value_3=0,rate_3=(-0.0005),lag_time_3=0.1){value_3=lag_time_3+rate_3*t
  if(value_3<min_value_3)value_3=min_value_3
  value_3
  },   
  
  #attack rate of adult predators on forage fish (used to be 0.1)
  d = 0.5       #Stocking rate for forage fish
)


#This is the starting conditions for the model. P is the predatory fish density,
#F is the forage fish density, and J is juvenile fish density
init_cond = c(P = 77, F = 0.067, J = 9.37)

#This function is the R code version of the system of equations that's described
#in the project description. It takes arguments t for time, y for current state
#of the community (defined here as the current densities of predators, forage
#fish and juvenile fish) and parms for parameters, that define how the
#population will change over time. It returns a rate of change for each of the
#state variables, so the ode function can use those rates of change to simulate
#how the population shifts over time
troph_tri_static_3 = function(t,y, parms){
  
  #This code extracts the three state variables (P,F, and J) from the y vector
  #so they can be referred to as single letters in the equations for rates of
  #change
  P = y["P"]
  F = y["F"]
  J = y["J"]
  
  #This next code calculates the derivatives at each point in time. 
  #the with(x,...) function here make the model parameters available by name
  #without having to type parms$e*P + parms$s...
  dP = with(parms, J/T_mat - m*P - e(t)*P + s) 
  dF = with(parms, r*F - b*F^2 - a_PF(t)*P*F + d)
  dJ = with(parms, f*P - J/T_mat - m*J - a_PJ*P*J - a_FJ*F*J)
  return(list(c(P=dP,F=dF, J=dJ)))
}

Time_plot_3 = seq(0,600,length.out = 1000)

#This code runs the model (it's an ordinary differential equation; thus, ode)
#the times argument says at what times you want to evaluate the model (here,
#from time 0 to time 200, split into 1000 equal steps). func is the function
#that calculates how quickly the population will change at each time point
#parms are the model parameters. see ?deSolve::ode for more information.
sim3 = deSolve::ode(y=init_cond,
                    times = Time_plot_3,
                    func = troph_tri_static_3,
                    parms = model_parameters_3)

plot(sim3)

time_series_e_3 <- sapply(Time_plot_3, model_parameters_3$e)
plot(time_series_e_3, type="l")

time_series_aPF_3 <- sapply(Time_plot_3, model_parameters_3$a_PF)
plot(time_series_aPF_3, type="l")



library(dplyr)

n_step_3 = length(Time_plot_3)

calc_1st_deriv_3 = function(y,delta_3) (lead(y,1) - lag(y,1))/(2*delta_3)
calc_2nd_deriv_3 = function(y,delta_3) (lead(y,1) + lag(y,1)-2*y)/delta_3^2

delta_3 <-sim3[2,"time"] - sim3[1,"time"]

first_deriv_3 <- matrix(NA,nrow= n_step_3, ncol =3)
first_deriv_3[,1] <- calc_1st_deriv_3(sim3[,"P"], delta_3)
first_deriv_3[,2] <- calc_1st_deriv_3(sim3[,"F"], delta_3)
first_deriv_3[,3] <- calc_1st_deriv_3(sim3[,"J"], delta_3)

second_deriv_3 <- matrix(NA,nrow= n_step_3, ncol =3)
second_deriv_3[,1] <- calc_2nd_deriv_3(sim3[,"P"], delta_3)
second_deriv_3[,2] <- calc_2nd_deriv_3(sim3[,"F"], delta_3)
second_deriv_3[,3] <- calc_2nd_deriv_3(sim3[,"J"], delta_3)

fisher_info_3 <- matrix(NA,nrow= n_step_3,ncol=1)

for(i in 1:n_step_3){
  numerator_3 <-  sum(first_deriv_3[i,]*second_deriv_3[i,])^2
  denominator_3 <- sqrt(sum(second_deriv_3[i,]^2))^6
  fisher_info_3[i,] <-  numerator_3/denominator_3 
}

library(tibbletime)
rolling_mean_3 <- rollify(mean, window = 100)
rolling_mean_fisher_3 <- rolling_mean_3(fisher_info_3[,1])

plot(sim3[,"time"],
     rolling_mean_fisher_3,
     type="l",
     ylab="Fisher Information - Rolling Mean",
     xlab="time", log='y')

troph_tri_jacobian_3 <- function(t,y,parms){
  
  jacobian_3 <- matrix(NA, nrow=3, ncol=3)
  
  #dP - differentiation variable: P 
  jacobian_3[1,1] <- with(model_parameters_3, (-e(t)-m))
  #dP - differentiation variable: F 
  jacobian_3[1,2] <- (0)
  #dP - differentiation variable: J 
  jacobian_3[1,3] <- with(model_parameters_3, (1/T_mat))
  
  
  #dF - differentiation variable: P 
  jacobian_3[2,1] <- with(model_parameters_3, (-a_PF(t)*y[2]))
  #dF - differentiation variable: F - error
  jacobian_3[2,2] <- with(model_parameters_3, ((-2*b*y[2])-a_PF(t)*y[1])+(r))
  #dF - differentiation variable: J
  jacobian_3[2,3] <- (0)
  
  
  #dJ - differentiation variable: P 
  jacobian_3[3,1] <- with(model_parameters_3, (f-(a_PJ*y[3])))
  #dJ - differentiation variable: F 
  jacobian_3[3,2] <- with(model_parameters_3, (-(a_FJ*y[3])))
  #dJ - differentiation variable: J 
  jacobian_3[3,3] <- with(model_parameters_3, ((-1/T_mat)-(a_PJ*y[1])-(m)-(a_FJ*y[2])))
  
  
  
  
  return(jacobian_3)
}

#create a vector of 600 time steps
times <- seq(0, 600, by=1)
n_steps <- length(times)

#create a data frame to hold equilibrium values
stable_states_3 <- data.frame(P = rep(0,times = n_steps),
                              F = rep(0,times = n_steps),
                              J = rep(0,times = n_steps),
                              eigen  = rep(0,times = n_steps),
                              time = times)

#set the current state for the loop to the original initial condition
current_state_3 <- init_cond

for(i in 1:n_steps){
  current_time <- times[i]
  #calculate the closest equilibrium point at the current time step
  root_value_3 <- stode(y= current_state_3,
                        time =current_time,
                        func = troph_tri_static_3,
                        jacfunc = troph_tri_jacobian_3,
                        parms = model_parameters_3,
                        positive = TRUE #this ensures that rootSolve will only find positive (or zero) solutions
  )
  
  #change the current state to this value
  current_state_3 <- root_value_3$y
  stable_states_3[i, c("P","F","J")] <- current_state_3
  
  #Calculate the Jacobian of the system at this equilibrium
  current_jacobian_3 <- troph_tri_jacobian_3(t = current_time, 
                                             y = current_state_3,
                                             parms = model_parameters_3)
  
  #calculate eigenvalues of this Jacobian and find the maximum real eigenvalue
  current_eigs_3 <- eigen(current_jacobian_3)
  stable_states_3[i,"eigen"] <- max(Re(current_eigs_3$values))
  
  #add a small perturbation to the current state to keep rootSolve from finding
  #only zero values after the regime shift.
  current_state_3 <- current_state_3 +0.5
}

#Find the regime shift point as the place where the eigen value of the Jacobian goes to zero (or just above)
regime_shift_3 <- stable_states_3$time[stable_states_3$eigen==max(stable_states_3$eigen)]

#time step of regime shift
regime_shift_3

#Plot all three time series, and the eigenvalues and Fisher information indices;
#mfcol specifies to add plots column-wise (fillin)
par(mfcol =c(3,2))

#plot the simulation, stable state, and annotate the regime shift for predators
plot(P~time, data = sim3, type="l")
points(P~time, data= stable_states_3,type="l",col ="red")
abline(v = regime_shift_3, col ="blue", lty=2)

#plot the simulation, stable state, and annotate the regime shift for forage fish
plot(F~time, data= sim3,type="l")
points(F~time, data= stable_states_3,type="l",col ="red")
abline(v = regime_shift_3, col ="blue", lty=2)

#plot the simulation, stable state, and annotate the regime shift for juveniles
plot(J~time, data= sim3,type="l")
points(J~time, data= stable_states_3,type="l",col ="red")
abline(v = regime_shift_3, col ="blue", lty=2)

#plotting eigenvalues
plot(eigen~time, data= stable_states_3,type="l")
abline(v = regime_shift_3, col ="blue", lty=2)

#plotting Fisher Information
plot(sim3[,"time"],
     rolling_mean_fisher_3,
     type="l",
     ylab="Fisher Information - Rolling Mean",
     xlab="time", log='y')
abline(v = regime_shift_3, col ="blue", lty=2)

par(mfcol =c(4,2))


