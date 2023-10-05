# Model files for chytrid model markdown document FitChytridModel.Rmd
# This file includes:
# 1. The ODE model
# 2. The parameter values and initial guesses (for beta, alpha, omega, and mu), 
#    and initial conditions
# 3. The function to calculate residuals using the ODE model function

# Trish Campbell & Jen Flegg & Claire Miller Sept 2023

# 1. ODE Model ---------------------------------------------------------------
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


# 2. Parameters --------------------------------------------------------------
# Rate parameters, given in Table S3
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
                mu = 0.0) # initial guess, to fit
# Initial conditions
x0 <- c( S=5,     I1=5,  I2=0,  I3=0, 
         SV=5, IV1=5, IV2=0, IV3=0)
x0["N"] = sum(x0)


# 3. Residuals function ------------------------------------------------------
# Calculates the residuals using the chytrid.ode.model function above
# inputs: 
#   - list of parameters (par) to fit with named elements called 
#     alpha, beta_greenhouse, beta_shaded, omega, and mu
#   - experimental data to fit to (expdata)
calc_residuals=function(par, expdata){
  calc_residuals_perenvironment <- function(par, Environment) {
    # times: match the experimental data
    t=unique(expdata$week) %>% sort()
    # Combine guessed parameters with known parameters
    # Add guessed parameters to variable containing parameters from above
    parameters['beta']=par[['beta']]
    parameters['alpha']=par[['alpha']]
    parameters['omega']=par[['omega']]
    parameters['mu']=par[['mu']]
    
    # solve ODE for a given set of parameters
    # initial concentrations from x0 (above)
    # function from above: chytrid.sir.model
    pred = ode(y=x0,times=t,func=chytrid.ode.model,parms=parameters) %>%
      data.frame()
    
    # Calculate residuals
    pred_long = pred %>% 
      pivot_longer(cols=-time,values_to="prediction")
    exp_long = filter(expdata,shaded==(Environment=="Shaded")) %>%
      select(-shaded,-Mesocosm) %>%
      pivot_longer(cols=-week,values_to="experiment") %>%
      rename(time="week")
    res = left_join(pred_long,exp_long,by=c("time","name")) %>%
      mutate(residuals=prediction-experiment)
    
    # Return residuals
    return(res$residuals)
  }
  
  # Run residuals function on greenhouse data
  par[["beta"]]=par[["beta_greenhouse"]]
  res_greenhouse = calc_residuals_perenvironment(par,"Greenhouse")
  # Run residuals function on shaded data
  par[["beta"]]=par[["beta_shaded"]]
  res_shaded = calc_residuals_perenvironment(par,"Shaded")
  
  # Return residuals from both environments
  return(c(res_greenhouse,res_shaded))
}
