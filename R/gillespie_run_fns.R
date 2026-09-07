
# 1. Function to run the Gillespie algorithm ---------------------------------
# @param parms: a named vector of parameters
# @param tf: total time
# @param Nsims: number of simulations
library(GillespieSSA)
run_gillespie <- function(parms, tf = 15, Nsims = 1)
{
  # Get the reactions and stochiometric matrix using function below
  gillespise_system <- build_stochastic_system()
  
  # Run the SSAs trajectories
  #-------------------------------------------------------------------
  # We don't explicitly model N so remove from initial conditions
  x0.gillespie <- x0[names(x0) != "N"]
  # Extract the final time
  
  # Run Nsims trajectories using ssa direct Gillespie method
  trajectories <- list()
  for (i in 1:Nsims) {
    trajectories[[i]] <- ssa(x0 = x0.gillespie,
                             a = gillespise_system[["a"]],
                             nu = gillespise_system[["nu"]],
                             parms = parms, tf = tf, 
                             method = ssa.d(),verbose = FALSE, consoleInterval = 1)
  }
  
  # Return the trajectories
  return(trajectories)
}


# 2. Function to extract specified times -------------------------------------
# (from gillespie simulation)
# sample_times = c(0,1,2,4,6,8,10,15)
# get_timepoints <- function(ssa_res, t_out = sample_times) {
get_timepoints <- function(ssa_res, t_out) {
  data = ssa_res$data
  rowi = sapply(t_out,FUN=(\(x) max(which(data[,"t"] < x),1)))
  out = data.frame(t=t_out, data[rowi, (colnames(data) != "t")])
  return(out)
}



# 3. Function to run gillespie for EasyABC ---------------------------------
run_gillespie_simulation_easyABC <- function(x_params, parm_names, base_params) {
    # Get the parameters and add the varying parameters
    parms_i <- base_params
    parms_i[parm_names] = x_params
    # Run the trajectory for shaded and unshaded
    parms_i["beta"] = parms_i["beta_sh"]
    res_shaded = run_gillespie(parms_i, Nsims=1)[[1]]
    parms_i["beta"] = parms_i["beta_un"]
    res_unshaded = run_gillespie(parms_i, Nsims=1)[[1]]
    # Return shaded and unshaded
    return(list("shaded"=res_shaded, "unshaded"=res_unshaded))
}