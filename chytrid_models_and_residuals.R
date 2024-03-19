# Model and residuals functions for chytrid model RMarkdown document FitChytridModel.Rmd
# This file includes:
# 1. The ODE model for a single environment
# 2. The ODE model for both environments combined
# 3. The parameter values and initial guesses (for beta, alpha, omega, and mu), 
#    and initial conditions
# 4. The function to calculate residuals for least squares
# 5. The function to calculate residuals for truncated normal observation model
# 6. The function to run the Gillespie algorithm for stochastic simulations

# Claire Miller & Jen Flegg  & Trish Campbell Feb 2024

# 1. ODE Model single environment --------------------------------------------
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

# 2. ODE Model both environments ---------------------------------------------
# inputs: 
#   - time (t, scalar), 
#   - current state vector (x, numeric vector), and 
#   - parameters (params, named numeric vector)
chytrid.model.envcomb <- function(t, x, params){
  p = params
  p["beta"] <- params["beta_un"]
  x.un <- x[grepl("_un",names(x))]
  y.un <- chytrid.ode.model(t, x.un, p)
  p["beta"] <- params["beta_sh"]
  x.sh <- x[grepl("_sh",names(x))]
  y.sh <- chytrid.ode.model(t, x.sh, p)
  return(list(c(y.un[[1]], y.sh[[1]])))
}

# 3. Parameters --------------------------------------------------------------
# Rate parameters, given in Table S3
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
x0["N"] = sum(x0)


# 3. Residuals function for least squares ------------------------------------
# Calculates the residuals using the chytrid.ode.model function above
# inputs: 
#   - list of parameters (par) to fit with named elements called 
#     alpha, beta_un, beta_sh, omega, and mu
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
  
  # Run residuals function on unshaded (greenhouse) data
  par[["beta"]]=par[["beta_un"]]
  res_greenhouse = calc_residuals_perenvironment(par,"Greenhouse")
  # Run residuals function on shaded data
  par[["beta"]]=par[["beta_sh"]]
  res_shaded = calc_residuals_perenvironment(par,"Shaded")
  
  # Return residuals from both environments
  return(c(res_greenhouse,res_shaded))
}

# 4. Residuals function for truncated normal ---------------------------------
# Calculates the residuals using the chytrid.ode.model function above
# inputs: 
#   - list of parameters (par) to fit with named elements called 
#     alpha, beta_greenhouse, beta_shaded, omega, and mu
#   - experimental data to fit to (expdata)
library(msm) # For the truncated normal distribution
calc_residuals_obsmodel=function(par, expdata){
  # Get time points from experimental data
  t=unique(expdata$week) %>% sort()
  
  # Add the unknowns to the parameters
  parameters['beta_un']=par[['beta_un']]
  parameters['beta_sh']=par[['beta_sh']]
  parameters['alpha']=par[['alpha']]
  parameters['omega']=par[['omega']]
  parameters['mu']=par[['mu']]
  
  # solve ODE for a given set of parameters
  x0.env = c(setNames(x0,paste0(names(x0),"_un")), 
             setNames(x0,paste0(names(x0),"_sh")))
  pred = ode(y=x0.env,times=t,
             func=chytrid.model.envcomb,
             parms=parameters) %>%
    data.frame()
  
  # Process experimental data
  pred_long = pred %>% filter(time>0) %>%
    pivot_longer(cols=-time,values_to="prediction")
  exp_long = expdata %>% filter(week>0) %>%
    select(-Mesocosm) %>%
    pivot_longer(cols=c(-week,-shaded),values_to="experiment") %>%
    rename(time="week") %>%
    mutate(name=ifelse(shaded,paste0(name,"_sh"),paste0(name,"_un"))) 
  
  # Calculate residuals
  res <- left_join(pred_long,exp_long,by=c("time","name"))
  negLL = -sum(dtnorm(res$experiment, res$prediction, 
                      sd = 1, lower = 0, log = TRUE))
  
  # Return residuals
  return(negLL)
}

# 5. Function to run the Gillespie algorithm ---------------------------------
# @param parms: a named vector of parameters
# @param tf: total time
# @param Nsims: number of simulations
library(GillespieSSA)
run_gillespie <- function(parms, tf = 15, Nsims = 10)
{
  # Build the reactions and stochiometric matrix (state-change matrix) 
  #-------------------------------------------------------------------
  # where each colums is a posible reaction with associate propensity in a vector
  # Matrix rows: S, I1, I2, I3, SV, IV1, IV2, IV3
  # we have 16 reactions and 8 states
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
  
  
  # Run the SSAs trajectories
  #-------------------------------------------------------------------
  # We don't explicitly model N so remove from initial conditions
  x0.gillespie <- x0[names(x0) != "N"]
  # Run Nsims trajectories using ssa direct Gillespie method
  trajectories <- list()
  for (i in 1:Nsims) {
    trajectories[[i]] <- ssa(x0 = x0.gillespie,
                             a = a,nu = nu,
                             parms = parms, tf = tf, 
                             method = ssa.d(),verbose = FALSE, consoleInterval = 1)
  }
  
  # Return the trajectories
  return(trajectories)
}