#Simulation 1 - Slowly increasing exploitation rate


library(deSolve)
library(ggplot2)
library(rootSolve)
library(mgcv)
library(dplyr)

#set model parameters

model_parameters_1 = list(
  T_mat = 5,    #length of time it takes a juvenile predator to mature to an adult
  m = 0.025,    #mortality rate of adult and juvenile fish
  s = 0.05,     #The stocking rate for adult predators
  e = function(t,min_value_1=0,max_value_1=1,rate_1=0.001,lag_time_1=0.15){value_1=rate_1*t-lag_time_1
  if(value_1<min_value_1)value_1=min_value_1
  if(value_1>max_value_1)value_1=max_value_1
  value_1
  },     
  
  #The extraction rate at any given time
  
  f = 0.5,      #amount of new offspring for each adult predator per unit time
  a_PJ  = 0.05, #Cannibalism rate of adult predators on juveniles
  a_FJ  = 0.1,  #attack rate of forage fish on juvenile predators
  
  r = 0.25,     #population growth rate of forage fish at low densities
  b = 0.005,    #density-dependence term for the forage fish
  a_PF = 0.1,   #attack rate of adult predators on forage fish
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
troph_tri_static_1 = function(t,y, parms){
  
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
  dF = with(parms, r*F - b*F^2 - a_PF*P*F + d)
  dJ = with(parms, f*P - J/T_mat - m*J - a_PJ*P*J - a_FJ*F*J)
  return(list(c(P=dP,F=dF, J=dJ)))
}

Time_plot_1 = seq(0,600,length.out = 100)

#This code runs the model (it's an ordinary differential equation; thus, ode)
#the times argument says at what times you want to evaluate the model (here,
#from time 0 to time 200, split into 1000 equal steps). func is the function
#that calculates how quickly the population will change at each time point
#parms are the model parameters. see ?deSolve::ode for more information.
sim1 = deSolve::ode(y=init_cond,
                    times = Time_plot_1,
                    func = troph_tri_static_1,
                    parms = model_parameters_1)


# Adding noise to simulations #### 

set.seed(10)
obs_error = 0.1

sim1_witherror = sim1

sim1_witherror[,"P"] = sim1[,"P"]*exp(rnorm(100,0,sd = obs_error))
sim1_witherror[,"J"] = sim1[,"J"]*exp(rnorm(100,0,sd = obs_error))
sim1_witherror[,"F"] = sim1[,"F"]*exp(rnorm(100,0,sd = obs_error))



plot(sim1_witherror)

sim1_P_gam = gam(P~s(time, k= 20),data= as.data.frame(sim1_witherror),
                 family = Gamma(link = "log"),method="REML")


plot(sim1_P_gam)

sim1_predicted = sim1_witherror

sim1_predicted[,"P"] = fitted(sim1_P_gam)


plot(sim1_predicted)


time_series_e_1 <- sapply(Time_plot_1, model_parameters_1$e)
plot(time_series_e_1, type="l")



n_step_1 = length(Time_plot_1)

calc_1st_deriv_1 = function(y,delta_1) (lead(y,1) - lag(y,1))/(2*delta_1)
calc_2nd_deriv_1 = function(y,delta_1) (lead(y,1) + lag(y,1)-2*y)/delta_1^2

delta_1 <-sim1[2,"time"] - sim1[1,"time"]

first_deriv_1 <- matrix(NA,nrow= n_step_1, ncol =3)
first_deriv_1[,1] <- calc_1st_deriv_1(sim1[,"P"], delta_1)
first_deriv_1[,2] <- calc_1st_deriv_1(sim1[,"F"], delta_1)
first_deriv_1[,3] <- calc_1st_deriv_1(sim1[,"J"], delta_1)

second_deriv_1 <- matrix(NA,nrow= n_step_1, ncol =3)
second_deriv_1[,1] <- calc_2nd_deriv_1(sim1[,"P"], delta_1)
second_deriv_1[,2] <- calc_2nd_deriv_1(sim1[,"F"], delta_1)
second_deriv_1[,3] <- calc_2nd_deriv_1(sim1[,"J"], delta_1)

fisher_info_1 <- matrix(NA,nrow= n_step_1,ncol=1)

for(i in 1:n_step_1){
  numerator_1 <-  sum(first_deriv_1[i,]*second_deriv_1[i,])^2
  denominator_1 <- sqrt(sum(second_deriv_1[i,]^2))^6
  fisher_info_1[i,] <-  numerator_1/denominator_1 
}

library(tibbletime)
rolling_mean_1 <- rollify(mean, window = 100)
rolling_mean_fisher_1 <- rolling_mean_1(fisher_info_1[,1])

plot(sim1[,"time"],
     rolling_mean_fisher_1,
     type="l",
     ylab="Fisher Information - Rolling Mean",
     xlab="time", log='y')

troph_tri_jacobian_1 <- function(t,y,parms){

  jacobian_1 <- matrix(NA, nrow=3, ncol=3)
  
  #dP - differentiation variable: P (-e(t)-m) 
  jacobian_1[1,1] <- with(parms, (-e(t)-m))
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

troph_tri_jacobian_1_sim <- troph_tri_jacobian_1(t = sim1[1000,"time"],sim1[1000,c("P","F","J")],parms = model_parameters_1)

eigen_jacobian_1 <- eigen(troph_tri_jacobian_1_sim)$values
eigen_jacobian_1
#eigen values
#eigen_jacobian_1$values
#eigen vectors
#eigen_jacobian_1$vectors

#create a vector of 600 time steps
times <- seq(0, 600, by=1)
n_steps <- length(times)

#create a data frame to hold equilibrium values
stable_states_1 <- data.frame(P = rep(0,times = n_steps),
                            F = rep(0,times = n_steps),
                            J = rep(0,times = n_steps),
                            eigen  = rep(0,times = n_steps),
                            time = times)

#set the current state for the loop to the original initial condition
current_state_1 <- init_cond

for(i in 1:n_steps){
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
regime_shift_1 <- stable_states_1$time[stable_states_1$eigen==max(stable_states_1$eigen)]

#time step of regime shift
regime_shift_1

#Plot all three time series, and the eigenvalues and Fisher information indices;
#mfcol specifies to add plots column-wise (fillin)
par(mfcol =c(3,2))

#plot the simulation, stable state, and annotate the regime shift for predators
plot(P~time, data = sim1, type="l")
points(P~time, data= stable_states_1,type="l",col ="red")
abline(v = regime_shift_1, col ="blue", lty=2)

#plot the simulation, stable state, and annotate the regime shift for forage fish
plot(F~time, data= sim1,type="l")
points(F~time, data= stable_states_1,type="l",col ="red")
abline(v = regime_shift_1, col ="blue", lty=2)

#plot the simulation, stable state, and annotate the regime shift for juveniles
plot(J~time, data= sim1,type="l")
points(J~time, data= stable_states_1,type="l",col ="red")
abline(v = regime_shift_1, col ="blue", lty=2)

#plotting eigenvalues
plot(eigen~time, data= stable_states_1,type="l")
abline(v = regime_shift_1, col ="blue", lty=2)

#plotting Fisher Information
plot(sim1[,"time"],
     rolling_mean_fisher_1,
     type="l",
     ylab="Fisher Information - Rolling Mean",
     xlab="time", log='y')
abline(v = regime_shift_1, col ="blue", lty=2)




