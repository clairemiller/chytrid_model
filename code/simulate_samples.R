
# Setup ------------------------------------------------------------------
dataset = "exp" # One of "syn" or "exp"

# Check command line arguments if they exist
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  dataset = args[1]
  check = assert_that(args %in% c("exp","syn"), 
                      msg="Arguments to script must be one of 'exp' or 'syn'.")
}
check = assert_that(length(args) <= 1, msg="Too many arguments passed to script.")

# Filepaths
output_data_filename <- paste0("../data/",dataset,"_priors_trajectories.RData")

# Generate the priors parameter sets --------------------------------------
set.seed(47)
priors = generate_priors(N_priors)
priors$iter_id <- 1:N_priors
mu = 0.021
if (dataset=="syn") {
  mu = lsfitparams[["mu"]]
}
priors$mu = mu

# Run trajectories ---------------------------------------------------------------
# Run in parallel
registerDoParallel(ncores)
tic("ABC simulate time")
simulation_list <- foreach(iter=1:nrow(priors)) %dopar% {
    # Get the prior for this iteration and merge with 'known'
    prior_i = priors[iter,]
    parms_i <- parameters 
    parms_i[names(prior_i)] <- prior_i 
    # Run the trajectory for shaded and unshaded
    parms_i["beta"] = parms_i["beta_sh"]
    res_shaded = run_gillespie(parms_i, Nsims=1)[[1]]
    parms_i["beta"] = parms_i["beta_un"]
    res_unshaded = run_gillespie(parms_i, Nsims=1)[[1]]
    # Get the correct time points for the simulated data
    sim_shaded = get_timepoints(res_shaded, sample_times)
    sim_unshaded = get_timepoints(res_unshaded, sample_times)
    # Add simulation details, combine, and return
    iter_id = prior_i$iter_id[1]
    sim_i = rbind(
      data.frame(sim_shaded, shaded=T, iter_id),
      data.frame(sim_unshaded, shaded=F, iter_id)
    )
    return(sim_i)
}
toc()


# Save data ---------------------------------------------------------------
save(simulation_list, priors, file=output_data_filename)
cat("Saved priors and trajectories to", output_data_filename, "\n")













