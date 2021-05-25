#New model from Boettiger and Batt

library(deSolve)
library(ggplot2)
library(rootSolve)

#Set the model parameters for the system

model_parameters_4 = list(
q = function(t, rate_4=(0.002/1020),init_q_4=0.001){
value_4 = t*rate_4 - init_q_4
return(value_4)
}, #harvest or stocking
s = 0.6, #survival rate fo juveniles to maturation
m_A = 0.4, #mortality rate of adults
f = 2, #fecundity of adult bass
c_FA = 0.3, #predation of planktivores by adult bass
c_JA = 0.1, #predation of juvenile bass by adult bass
c_JF = 0.5, #predation of juvenile bass by planktivorous fish
F_o = 200, #abundance of planktivorous fish in non-foraging arena
D_F = 0.09, #diffusion of planktivores between refuge and forage arena
v = 80, #rate at which J enter a foraging arena and become vulnerable
h = 80 #rate at which J hide in a refuge
)

#adjusted model parameter values from the paper
#model_parameters_4 = list(
  #q = function(t, rate_4=0.00000196078,init_q_4=0.001){
    #value_4 = t*rate_4 - init_q_4
 #return(value_4)
#}, #harvest or stocking
  #s = 0.5, #survival rate fo juveniles to maturation
  #m_A = 0.1, #mortality rate of adults
  #f = 2, #fecundity of adult bass
  #c_FA = 0.3, #predation of planktivores by adult bass
  #c_JA = 0.001, #predation of juvenile bass by adult bass
  #c_JF = 0.5, #predation of juvenile bass by planktivorous fish
  #F_o = 100, #abundance of planktivorous fish in non-foraging arena
  #D_F = 0.1, #diffusion of planktivores between refuge and forage arena
  #v = 1, #rate at which J enter a foraging arena and become vulnerable
  #h = 8 #rate at which J hide in a refuge
#)

init_cond_4 = c(A = 12, F = 1, J = 8)

#set the function for the system equations

trophic_system_4 = function(t,y, parms){
  
  A = y["A"]
  F = y["F"]
  J = y["J"]
  
  dA = with(parms, s*J - q(t)*A - m_A*A)
  dF = with(parms, D_F*(F_o - F) - c_FA*F*A)
  dJ = with(parms, f*A - c_JA*J*A - (c_JF*v*J*F)/(h + v + c_JF*F) - s*J)
  return(list(c(A=dA,F=dF, J=dJ)))
} 

#set time for plot

Time_plot_4 = seq(0,1020,length.out = 102000)

#plot the simulation

sim4 = deSolve::ode(y=init_cond_4,
                    times = Time_plot_4,
                    func = trophic_system_4,
                    parms = model_parameters_4)

plot(sim4)

time_series_q_4 <- sapply(Time_plot_4, model_parameters_4$q)
plot(time_series_q_4, type="l")


#model_parameters_4_1 = list(
  #q = function(t, rate_4_1=0.001,init_q_4_1=0.001){
    #value_4_1 = t*rate_4_1 - init_q_4_1
    #return(value_4_1)
  #}, #harvest or stocking
  #s = 0.6, #survival rate fo juveniles to maturation
  #m_A = 0.4, #mortality rate of adults
  #f = 2, #fecundity of adult bass
  #c_FA = 0.3, #predation of planktivores by adult bass
  #c_JA = 0.1, #predation of juvenile bass by adult bass
  #c_JF = 0.5, #predation of juvenile bass by planktivorous fish
  #F_o = 200, #abundance of planktivorous fish in non-foraging arena
  #D_F = 0.09, #diffusion of planktivores between refuge and forage arena
  #v = 80, #rate at which J enter a foraging arena and become vulnerable
  #h = 80 #rate at which J hide in a refuge
#)

#sim4_1 = deSolve::ode(y=init_cond_4,
                    #times = Time_plot_4,
                    #func = trophic_system_4,
                    #parms = model_parameters_4_1)

#plot(sim4_1)

#time_series_q_4_1 <- sapply(Time_plot_4, model_parameters_4_1$q)
#plot(time_series_q_4_1, type="l")


library(dplyr)

#set up fisher informaiton calculations for the simulation

n_step_4 = length(Time_plot_4)

calc_1st_deriv_4 = function(y,delta_4) (lead(y,1) - lag(y,1))/(2*delta_4)
calc_2nd_deriv_4 = function(y,delta_4) (lead(y,1) + lag(y,1)-2*y)/delta_4^2

delta_4 <-sim4[2,"time"] - sim4[1,"time"]

first_deriv_4 <- matrix(NA,nrow= n_step_4, ncol =3)
first_deriv_4[,1] <- calc_1st_deriv_4(sim4[,"A"], delta_4)
first_deriv_4[,2] <- calc_1st_deriv_4(sim4[,"F"], delta_4)
first_deriv_4[,3] <- calc_1st_deriv_4(sim4[,"J"], delta_4)

second_deriv_4 <- matrix(NA,nrow= n_step_4, ncol =3)
second_deriv_4[,1] <- calc_2nd_deriv_4(sim4[,"A"], delta_4)
second_deriv_4[,2] <- calc_2nd_deriv_4(sim4[,"F"], delta_4)
second_deriv_4[,3] <- calc_2nd_deriv_4(sim4[,"J"], delta_4)

fisher_info_4 <- matrix(NA,nrow= n_step_4,ncol=1)

for(i in 1:n_step_4){
  numerator_4 <-  sum(first_deriv_4[i,]*second_deriv_4[i,])^2
  denominator_4 <- sqrt(sum(second_deriv_4[i,]^2))^6
  fisher_info_4[i,] <-  numerator_4/denominator_4 
}

library(tibbletime)
rolling_mean_4 <- rollify(mean, window = 100)
rolling_mean_fisher_4 <- rolling_mean_4(fisher_info_4[,1])

plot(sim4[,"time"],
     rolling_mean_fisher_4,
     type="l",
     ylab="Fisher Information - Rolling Mean",
     xlab="time", log='y')

#set up jacobian matrix

troph_tri_jacobian_4 <- function(t,y,parms){
  
  jacobian_4 <- matrix(NA, nrow=3, ncol=3)
  
  #dA - differentiation variable: A -q(t) - m_a 
  jacobian_4[1,1] <- with(parms, - q(t) - m_A)
  #dA - differentiation variable: F 
  jacobian_4[1,2] <- (0)
  #dA - differentiation variable: J 
  jacobian_4[1,3] <- with(parms, s)
  
  
  #dF - differentiation variable: A 
  jacobian_4[2,1] <- with(parms, -c_FA*y[2])
  #dF - differentiation variable: F 
  jacobian_4[2,2] <- with(parms, -D_F - y[1]*c_FA)
  #dF - differentiation variable: J 
  jacobian_4[2,3] <- (0)
  
  #dJ - differentiation variable: A 
  jacobian_4[3,1] <- with(parms, y[2] - c_JA*y[3])
  #dJ - differentiation variable: F
  jacobian_4[3,2] <- with(parms, -((c_JF*y[3]*v)/(c_JF*y[2]+v+h))+((((c_JF)^2)*y[3]*v*y[2])/(((c_JF*y[2])+v+h)^2))+y[1])
  #dJ - differentiation variable: J
  jacobian_4[3,3] <- with(parms, -((c_JF*y[2]*v)/(v+h+c_JF*y[2]))-s-(y[1]*c_JA))
  
  
  
  return(jacobian_4)
}

#troph_tri_jacobian_4_sim <- troph_tri_jacobian_4(t = sim4[1000,"time"],sim1[1000,c("A","F","J")],parms = model_parameters_4)

#eigen_jacobian_4 <- eigen(troph_tri_jacobian_1_sim)$values
#eigen_jacobian_4
#eigen values
#eigen_jacobian_1$values
#eigen vectors
#eigen_jacobian_1$vectors

#create a vector of 600 time steps
times_4 <- seq(0, 1020, by=1)
n_steps_4 <- length(times_4)

#create a data frame to hold equilibrium values
stable_states_4 <- data.frame(A = rep(0,times_4 = n_steps_4),
                              F = rep(0,times_4 = n_steps_4),
                              J = rep(0,times_4 = n_steps_4),
                              eigen  = rep(0,times_4 = n_steps_4),
                              time = times_4)

#set the current state for the loop to the original initial condition
current_state_4 <- init_cond_4

for(i in 1:n_steps_4){
  current_time_4 <- times_4[i]
  #calculate the closest equilibrium point at the current time step
  root_value_4 <- stode(y= current_state_4,
                        time =current_time_4,
                        func = trophic_system_4,
                        jacfunc = troph_tri_jacobian_4,
                        parms = model_parameters_4,
                        positive = TRUE #this ensures that rootSolve will only find positive (or zero) solutions
  )
  
  #change the current state to this value
  current_state_4 <- root_value_4$y
  stable_states_4[i, c("A","F","J")] <- current_state_4
  
  #Calculate the Jacobian of the system at this equilibrium
  current_jacobian_4 <- troph_tri_jacobian_4(t = current_time_4, 
                                             y = current_state_4,
                                             parms = model_parameters_4)
  
  #calculate eigenvalues of this Jacobian and find the maximum real eigenvalue
  current_eigs_4 <- eigen(current_jacobian_4)
  stable_states_4[i,"eigen"] <- max(Re(current_eigs_4$values))
  
  #add a small perturbation to the current state to keep rootSolve from finding
  #only zero values after the regime shift.
  current_state_4 <- current_state_4 +0.5
}

#Find the regime shift point as the place where the eigen value of the Jacobian goes to zero (or just above)
regime_shift_4 <- stable_states_4$time[stable_states_4$eigen==max(stable_states_4$eigen)]

#time step of regime shift
regime_shift_4

#Plot all three time series, and the eigenvalues and Fisher information indices;
#mfcol specifies to add plots column-wise (fillin)
par(mfcol =c(3,2))

#plot the simulation, stable state, and annotate the regime shift for predators
plot(A~time, data = sim4, type="l")
points(A~time, data= stable_states_4,type="l",col ="red")
abline(v = regime_shift_4, col ="blue", lty=2)

#plot the simulation, stable state, and annotate the regime shift for forage fish
plot(F~time, data= sim4,type="l")
points(F~time, data= stable_states_4,type="l",col ="red")
abline(v = regime_shift_4, col ="blue", lty=2)

#plot the simulation, stable state, and annotate the regime shift for juveniles
plot(J~time, data= sim4,type="l")
points(J~time, data= stable_states_4,type="l",col ="red")
abline(v = regime_shift_4, col ="blue", lty=2)

#plotting eigenvalues
plot(eigen~time, data= stable_states_4,type="l")
abline(v = regime_shift_4, col ="blue", lty=2)

#plotting Fisher Information
plot(sim4[,"time"],
     rolling_mean_fisher_4,
     type="l",
     ylab="Fisher Information - Rolling Mean",
     xlab="time", log='y')
abline(v = regime_shift_4, col ="blue", lty=2)



