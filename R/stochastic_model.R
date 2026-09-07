# Model and residuals functions for scripts: script_abc_rejection.R
# This file includes:
# 1. The model parameters and initial conditions (including ODE fit parameters)
# 2. Function that returns the stochastic system
#    (reactions and stochiometric matrix (state-change matrix))

# Claire Miller & Jen Flegg  & Trish Campbell Feb 2024

# 1. Parameters --------------------------------------------------------------
parameters <- c(beta = 0.006, # initial guess, to fit
                alpha = 0.3, # initial guess, to fit
                omega = 1/120, # initial guess, to fit
                gamma1 = 1.0/2.5, 
                gamma2 = 1.0/4.5, 
                h_int = 100,
                m_int = 10,
                h_intv = 10,
                m_intv = 1,
                l_intv = 0.1,
                mu = NA) # loss rate, estimated from data, updated in run_seq_abc.R for each setup

# Initial conditions
x0 <- c( S=5,  I1=5,  I2=0,  I3=0, 
         SV=5, IV1=5, IV2=0, IV3=0)

# Parameters determined using the LS method
lsfitparams <- c(
  beta_sh = 0.05,
  beta_un = 0.02,
  alpha = 0.177,
  omega = 0.258,
  mu = 0.029
)

# 2. Build the stochastic system ---------------------------------
# Builds the reactions and stochiometric matrix (state-change matrix) 
# where each colums is a posible reaction with associate propensity in a vector
# Matrix rows: S, I1, I2, I3, SV, IV1, IV2, IV3
# we have 16 reactions and 8 states
build_stochastic_system <- function() {

  a  <- rep(NA,16) # reactions
  nu <- matrix(0,  nrow=8,ncol = 16,byrow=TRUE) # stochiometric matrix
  # S->I
  a[1] <- "beta*S*(m_int * I1 + m_intv * IV1 + h_int * I2 + h_intv * IV2 + I3 + l_intv * IV3)/(S+I1+I2+I3+SV+IV1+IV2+IV3-1)"
  nu[c(1,2),1] <- c(-1,1)
  # S loss
  a[2] <- "mu*S"
  nu[c(1),2] <- c(-1)
  # I1->I2
  a[3] <- "gamma1*I1"
  nu[c(2,3),3] <- c(-1,1)
  # I1 loss
  a[4] <- "mu*I1"
  nu[c(2),4] <- c(-1)
  # I2->I3
  a[5] <- "gamma2*I2"
  nu[c(3,4),5] <- c(-1,1)
  # I2 loss
  a[6] <- "mu*I2"
  nu[c(3),6] <- c(-1)
  # I3->SV
  a[7] <- "omega*I3"
  nu[c(4,5),7] <- c(-1,1)
  # I3 loss
  a[8] <- "mu*I3"
  nu[c(4),8] <- c(-1)
  # SV->IV1
  a[9] <- "alpha*beta*SV*(m_int * I1 + m_intv * IV1 + h_int * I2 + h_intv * IV2 + I3 + l_intv * IV3)/(S+I1+I2+I3+SV+IV1+IV2+IV3-1)"
  nu[c(5,6),9] <- c(-1,1)
  # SV loss
  a[10] <- "mu*SV"
  nu[c(5),10] <- c(-1)
  # IV1->IV2
  a[11] <- "gamma1*IV1"
  nu[c(6,7),11] <- c(-1,1)
  # IV1 loss
  a[12] <- "mu*IV1"
  nu[c(6),12] <- c(-1)
  # IV2->IV3
  a[13] <- "gamma2*IV2"
  nu[c(7,8),13] <- c(-1,1)
  # IV2 loss
  a[14] <- "mu*IV2"
  nu[c(7),14] <- c(-1)
  # IV3->SV
  a[15] <- "omega*IV3"
  nu[c(5,8),15] <- c(1,-1)
  # IV3 loss
  a[16] <- "mu*IV3"
  nu[c(8),16] <- c(-1)
  
  # Return as list
  return(list("a"=a,"nu"=nu))
}

