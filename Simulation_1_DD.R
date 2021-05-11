#Simulation 1 - Slowly increasing exploitation rate


library(deSolve)
library(ggplot2)

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

Time_plot_1 = seq(0,600,length.out = 1000)

#This code runs the model (it's an ordinary differential equation; thus, ode)
#the times argument says at what times you want to evaluate the model (here,
#from time 0 to time 200, split into 1000 equal steps). func is the function
#that calculates how quickly the population will change at each time point
#parms are the model parameters. see ?deSolve::ode for more information.
sim1 = deSolve::ode(y=init_cond,
                    times = Time_plot_1,
                    func = troph_tri_static_1,
                    parms = model_parameters_1)

plot(sim1)

plot(model_parameters_1$e(t),t=0,600)

library(dplyr)

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
jacobian_1[1,1] <- with(model_parameters_1, (-e(t)-m))
#dP - differentiation variable: J (1/T_mat)
jacobian_1[2,1] <- with(model_parameters_1, (1/T_mat))
#dP - differentiation variable: F 0
jacobian_1[3,1] <- (0)

#dJ - differentiation variable: P (F-(a_PJ*J)
jacobian_1[1,2] <- with(model_parameters_1, (init_cond[2]-(a_PJ*init_cond[3])))
#dJ - differentiation variable: J (-1/T_mat)-(a_PJ*P)-(m)-(a_FJ*F)
jacobian_1[2,2] <- with(model_parameters_1, ((-1/T_mat)-(a_PJ*init_cond[1])-(m)-(a_FJ*init_cond[2])))
#dJ - differentiation variable: F (P-(a_FJ*J)
jacobian_1[3,2] <- with(model_parameters_1, (init_cond[1]-(a_FJ*init_cond[3])))

#dF - differentiation variable: P (-a_PF*F)
jacobian_1[1,3] <- with(model_parameters_1, (-a_PF*init_cond[2]))
#dF - differentiation variable: J 0
jacobian_1[2,3] <- (0)
#dF - differentiation variable: F (-2*b*F)+(r)-(a_PF*P)
jacobian_1[3,3] <- with(model_parameters_1, ((-2*b*init_cond[2])+(r)-(a_PF*init_cond[1])))

return(jacobian_1)
}

troph_tri_jacobian_1_sim <- troph_tri_jacobian_1(t = sim1[1000,"time"],sim1[1000,c("P","F","J")],parms = model_parameters_1)

eigen_jacobian_1 <- eigen(troph_tri_jacobian_1_sim)$values
eigen_jacobian_1
#eigen values
#eigen_jacobian_1$values
#eigen vectors
#eigen_jacobian_1$vectors




