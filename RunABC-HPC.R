# Command to run using singularity: 
#   singularity exec ../Containers/chytrid_R_container.sif Rscript RunABC.R

# Load libraries
library(magrittr) 
library(GillespieSSA) # To run gillespie trajectories
library(dplyr)
library(tidyr)
library(knitr)
library(msm) # For truncated normal distribution
# Parallel implementation
library(foreach)
library(doParallel) 

# Setup
#----------------------------------------------------------------------------------------------
iterations = 5e6
nprocs = 32
set.seed(1) # Set the seed (for gillespie)

# Read in data
obs.long <- read.csv("data/expdata.csv")
# Change the time column from week to t to match function nomenclature
obs <- rename(obs.long,"t"="week") %>%
# We want to match the wide format with one time column 't' and one column per compartment
    pivot_wider(names_from="compartment", values_from="N", values_fill = 0) %>%
    select(-Mesocosm)


# Functions to run the system
#----------------------------------------------------------------------------------------------
# Source the ODE system from the R script
source("chytrid_models_and_residuals.R")
# Get the gillespie system using the function from this script
gillespie_system <- build_system_gillespie()

# A function that runs the model for a certain set of params
# @param params the parameter guess we want to run
# @Tvec a vector of times we want to get the output for
runModel <- function(params.fit, Tvec)
{
  # Initial conditions
  x0 <- c( S=5,  I1=5,  I2=0,  I3=0, 
           SV=5, IV1=5, IV2=0, IV3=0)
  # Combine fitting parameters with other model parameters 
  # (These are defined in chytrid_models_and_residuals.R)
  parameters[names(params.fit)] = params.fit
  # Run a trajectory
  trajectory <- ssa(x0 = x0,
      a = gillespie_system[["a"]],
      nu = gillespie_system[["nu"]],
      parms = parameters, 
      tf = max(Tvec),
      method = ssa.d(),verbose = FALSE, consoleInterval = 1)
  # Extract the data at the correct time points
  Xi <- lapply(setNames(Tvec,Tvec),function(ti, data) {
    # Determine the appropriate row in the solution (max is in case ti=0)
    rowi = max(which(data[,"t"] >= ti)[1]-1,1)
    # If it returns NA it will be the last time point
    rowi = ifelse(is.na(rowi), nrow(data), rowi)
    return(data[rowi,])
  }, data=trajectory$data)
  Xi <- bind_rows(Xi, .id="t") # override the t column with obs time
  return(Xi)
}

runModelBothEnv <- function(params.fit, Tvec)
{
  # Environment independent parameters
  params.env <- params.fit[grep("beta",names(params.fit),invert=T)]
  # Run shaded model
  params.shaded = c(params.env, beta=params.fit[["beta_sh"]])
  Xsim_sh <- runModel(params.shaded, Tvec)
  # Run unshaded model
  params.unshaded = c(params.env, beta=params.fit[["beta_un"]])
  Xsim_un <- runModel(params.unshaded, Tvec)
  # Return both as dataframe with boolean column for shaded or not
  Xsim_sh$shaded = T
  Xsim_un$shaded = F
  return(rbind(Xsim_sh, Xsim_un))
}


# Functions to calculate the summary statistics
#----------------------------------------------------------------------------------------------
# Function to calculate the summary statistic for acceptance
# @param Xobs a dataframe of observed data with a time column named 't' and one column per compartment
# @param Xsim a dataframe of simulated data matching structure of Xobs
# @return S a named vector of summary statistics with one element per compartment
calcSummaryStatisticsOneEnv <- function(Xobs, Xsim) {
  # Check the dataframes are matching structures (colnames and time)
  stopifnot(
    colnames(Xsim) %in% colnames(Xobs), # Order apathetic
    colnames(Xobs) %in% colnames(Xsim), # Order apathetic
    unique(Xobs$t) == unique(Xsim$t)
  )
  # Now we need to calculate the differences for each time point
  # Split the observations into individual time points
  Xdiff <- split(Xobs, Xobs$t) %>%
     lapply(function(xobsi) {
       # Get the relevant simulation row and exclude time column
       xsimi <- Xsim[Xsim$t==xobsi$t[1],names(Xsim) != "t"] 
       # Calculate the difference between obs and sim
       # Note mapply is much faster than data frame subtraction
       # (xobsi[,names(xsimi)] ensures the same column order)
       sdiff <- mapply("-", xobsi[,names(xsimi)], xsimi)
       # Return the absolute values of the differences
       return((sdiff)^2)
     }) %>% 
     # Convert back to a data frame
     do.call(rbind, .)
  # Return the compartment sums / total number of observations
  return(sum(Xdiff))
}

calcSummaryStatisticsBothEnv <- function(Xobs, Xsim) {
  # Shaded
  Xobs_sh <- Xobs[Xobs$shaded, colnames(Xobs) != "shaded"]
  Xsim_sh <- Xsim[Xsim$shaded, colnames(Xsim) != "shaded"]
  Si_sh <- calcSummaryStatisticsOneEnv(Xobs_sh, Xsim_sh)
  # Unshaded
  Xobs_un <- Xobs[!Xobs$shaded, colnames(Xobs) != "shaded"]
  Xsim_un <- Xsim[!Xsim$shaded, colnames(Xsim) != "shaded"]
  Si_un <- calcSummaryStatisticsOneEnv(Xobs_un, Xsim_un)
  # Return the average of shaded and unshaded
  return((Si_sh+Si_un)/2)
}


# Define the priors and run the ABC (save results)
#----------------------------------------------------------------------------------------------
# First let's generate all the proposed parameter sets as a list
parametersets <- lapply(1:iterations, function(i) {
    c(beta_sh = runif(1, min = 0, max = 1),
      beta_un = runif(1, min = 0, max = 1),
      alpha = rtnorm(1, mean = 1, sd = 1, lower = 0),
      omega = runif(1, min = 0, max = 7),
      mu = runif(1, min = 0, max = 1) )
  })

# Get the vector of output times from observed data
tOut <- unique(obs$t)
# Run model in parallel for each parameter set
registerDoParallel(nprocs)
run.time <- system.time({
  XsimList <- foreach(p=parametersets) %dopar% {
    runModelBothEnv(params.fit = p, Tvec = tOut)
  }
})
run.time <- data.frame(as.list(run.time)) 
kable(run.time, caption=paste("Model run computational time for",iterations,"iterations (on",nprocs,"CPUs)."))
# Calculate the summary statistics
summ.time <- system.time({
SList <- foreach(Xsimi = XsimList) %dopar% {
       calcSummaryStatisticsBothEnv(obs, Xsimi)
}
})
summ.time <- data.frame(as.list(summ.time))
kable(summ.time, caption=paste("Summary statistic computational time for",iterations,"iterations (on",nprocs,"CPUs)."))

# Save the results (separately cos XsimList is quite big and we only expect to use SList)
save(parametersets, XsimList, file="data/ABCSimulationData.RData")
save(parametersets, SList, file="data/ABCSummaryStatistic.RData")
# Also output the computational time
sink("data/ABCompTime.txt")
kables(list(kable(run.time,caption="Model run"), kable(summ.time,caption="Summary statistic calculation")), 
        caption=paste("Computational time for",iterations,"iterations (on",nprocs,"CPUs)."))
sink()