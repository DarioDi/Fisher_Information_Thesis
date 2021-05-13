#New model from Boettiger and Batt

library(deSolve)
library(ggplot2)
library(rootSolve)

#Set the model parameters for the system

model_parameters_4 = list(
  q = function(t, rate_4=0.001,init_q_4=0.001){
    value_4 = t*rate_4 - init_q_4
 return(value_4)
}, #harvest or stocking
  s = 0.5, #survival rate fo juveniles to maturation
  m_A = 0.1, #mortality rate of adults
  f = 2, #fecundity of adult bass
  c_FA = 0.3, #predation of planktivores by adult bass
  c_JA = 0.001, #predation of juvenile bass by adult bass
  c_JF = 0.5, #predation of juvenile bass by planktivorous fish
  F_o = 100, #abundance of planktivorous fish in non-foraging arena
  D_F = 0.1, #diffusion of planktivores between refuge and forage arena
  v = 1, #rate at which J enter a foraging arena and become vulnerable
  h = 8 #rate at which J hide in a refuge
)

init_cond_4 = c(A = 1.5, F = 25, J = 2)

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


model_parameters_4_1 = list(
  q = function(t, rate_4_1=0.001,init_q_4_1=0.001){
    value_4_1 = t*rate_4_1 - init_q_4_1
    return(value_4_1)
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

sim4_1 = deSolve::ode(y=init_cond_4,
                    times = Time_plot_4,
                    func = trophic_system_4,
                    parms = model_parameters_4_1)

plot(sim4_1)

time_series_q_4_1 <- sapply(Time_plot_4, model_parameters_4_1$q)
plot(time_series_q_4_1, type="l")




