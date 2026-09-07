# Model and residuals functions for scripts:
# ?
# This file includes:
# 1. The model parameters and initial conditions
# 2. The ODE model for a single environment

# Claire Miller & Jen Flegg  & Trish Campbell Feb 2024

# Parameters --------------------------------------------------------------
parameters <- c(beta = 0.006, # initial guess, to fit
                alpha = 0.3, # initial guess, to fit
                omega = 1/120,
                gamma1 = 1.0/2.5, 
                gamma2 = 1.0/4.5, 
                h_int = 100,
                m_int = 10,
                h_intv = 10,
                m_intv = 1,
                l_intv = 0.1,
                mu = 0.0) # initial guess, to fit

# Initial conditions
x0 <- c( S=5,  I1=5,  I2=0,  I3=0, 
         SV=5, IV1=5, IV2=0, IV3=0)
# x0["N"] = sum(x0)

# Parameters determined using the LS method
lsfitparams <- c(
  beta_sh = 0.05,
  beta_un = 0.02,
  alpha = 0.177,
  omega = 0.258,
  mu = 0.029
)


# ODE Model  --------------------------------------------
# inputs: 
#   - time (t, scalar), 
#   - current state vector (x, numeric vector), and 
#   - parameters (params, named numeric vector)
chytrid.ode.model <- function (t, x, params) {
  # extract the state variables
  S   <- x[1]
  I1  <- x[2]
  I2  <- x[3]
  I3  <- x[4]
  SV  <- x[5]
  IV1 <- x[6]
  IV2 <- x[7]
  IV3 <- x[8]
  N   <- x[9]
  
  # extract the parameters
  beta   <- params["beta"] # transmission rate coefficient
  alpha  <- params["alpha"] # rate of infection in vaccinated relative to naive
  omega  <- params["omega"] # transition rate I3 to SV/IV3 to SV
  gamma1 <- params["gamma1"] # reciprocal of time spent in I1
  gamma2 <- params["gamma2"] # reciprocal of time spent in I2
  h_int  <- params["h_int"] # intensity in naive high group relative to low naive
  m_int  <- params["m_int"] # intensity in naive medium group relative to low naive
  h_intv <- params["h_intv"] # intensity in vacc high group relative to low naive
  m_intv <- params["m_intv"] # intensity in vacc medium group relative to low naive
  l_intv <- params["l_intv"] # intensity in vacc low group relative to low naive
  mu <- params["mu"] # mortality rate
  
  # weighted sum of infectious
  I <- m_int * I1 + m_intv * IV1 + h_int * I2 + h_intv * IV2 + I3 + l_intv * IV3  
  
  # the model equations
  dSdt   <- -beta*S*I/(N-1) - mu*S
  dI1dt  <- beta*S*I/(N-1) - (gamma1+mu)*I1
  dI2dt  <- gamma1*I1 - (gamma2+mu)*I2
  dI3dt  <- gamma2*I2 - (omega+mu)*I3
  dSVdt  <- -alpha*beta*SV*I/(N-1) + omega*IV3 + omega*I3 - mu*SV
  dIV1dt <- alpha*beta*SV*I/(N-1) - (gamma1+mu)*IV1
  dIV2dt <- gamma1*IV1 - (gamma2+mu)*IV2
  dIV3dt <- gamma2*IV2 - (omega+mu)*IV3
  dNdt <- -mu*(S+I1+I2+I3+SV+IV1+IV2+IV3)
  
  # combine results into a single vector
  dxdt <- c(dSdt, dI1dt, dI2dt, dI3dt, dSVdt, dIV1dt, dIV2dt, dIV3dt, dNdt)
  
  # return result as a list
  list(dxdt)
}