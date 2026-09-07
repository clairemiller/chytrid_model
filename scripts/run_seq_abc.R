# Based on the vignette: https://cran.r-project.org/web/packages/EasyABC/vignettes/EasyABC.html
library(EasyABC)
library(tidyr)
library(dplyr)
devtools::load_all(".")


# 1. Choose our dataset, algorithm hyperparameters, and summary statistic --------------
# Choose which study to run: synthetic (TRUE) or experimental (FALSE)
synthetic_study = TRUE
base_params = parameters
if (synthetic_study) {
  obs_data = syndata
  root_output_dir = "data/seq-abc_syn-results"
  base_params["mu"] = 0.033
} else {
  obs_data = expdata
  root_output_dir = "data/seq-abc_exp-results"
  base_params["mu"] = 0.021
}
print(base_params)

# Algorithm hyperparameters
abc_params = list("method" = "Lenormand", "alpha" = 0.4, 
                  "p_acc_min" = 0.005, "n_particles" =  5e4) # Actual, ~6-7 days
                  # "p_acc_min" = 0.5, "n_particles" =  100) # For testing ~5 seconds

# Summary statistic - this is selected in R/seq_abc_summ_stats.R
cat("Summary statistic function:", summary_stat_fn_name, "\n")

# Generate the output directory name based on the algorithm parameters
output_dir_name = paste0("easyABC_output-", summary_stat_fn_name, 
                        "-alpha", abc_params[["alpha"]], 
                        "-n_particles", abc_params[["n_particles"]])
output_path = file.path(root_output_dir, output_dir_name)
cat("Output path:", output_path, "\n")



# 2. Define the priors ---------------------------------------------------
priors <- list(
  "beta_un" = c("unif", 0, 1.0), # beta_un
  "beta_sh" = c("unif", 0, 1.0), # beta_sh
  "alpha" = c("unif", 0, 2.0), # alpha
  "omega" = c("unif", 0, 10.0) # omega
)



# 3. Define the simulation code --------------------------------------------
run_model <- function(x){
  parm_names = names(priors)
  stopifnot(length(x) == length(parm_names)) # This does not set a seed
  ssa_results = run_gillespie_simulation_easyABC(x, parm_names, base_params)
  sim_data <- match_simulations_to_observations(
                res_shaded = ssa_results[["shaded"]], 
                res_unshaded = ssa_results[["unshaded"]], 
                obs_data = obs_data) 
  summ_stat = summary_stat_fn(sim_data_df = sim_data, obs_data_df = obs_data)
  return(summ_stat[[1]]) # just return the summary statistic, not the labels
}

make_model_parallel <- function(obs_data, parm_names, base_params) {
  # To make the following data accessible to the parallel workers
  force(obs_data)
  force(parm_names)
  force(base_params)
  stopifnot(!is.na(base_params["mu"])) # Make sure the loss rate is set

  function(x) {
    devtools::load_all(".") # Load the model and functions for each worker
    set.seed(x[1]) # This sets a seed for reproducibility
    stopifnot(length(x) == length(parm_names) + 1) # to make sure this is using a seed
    ssa_results = run_gillespie_simulation_easyABC(x[-1], parm_names, base_params)
    # Calculate the summary statistic and return
    sim_data <- match_simulations_to_observations(
                        res_shaded = ssa_results[["shaded"]], 
                        res_unshaded = ssa_results[["unshaded"]], 
                        obs_data = obs_data) 
    summ_stat = summary_stat_fn(sim_data_df = sim_data, obs_data_df = obs_data) 
    return(summ_stat[[1]]) # just return the summary statistic, not the labels
  }
}
run_model_parallel <- make_model_parallel(obs_data, names(priors), base_params)



# 4. Calculate the summary statistics for the observed data ---------------
summ_stat_obs = summary_stat_fn(obs_data_df = obs_data, sim_data_df = obs_data)



# 5. Run ABC Rejection ------------------------------------------------------
project_dir = getwd() # Save the current working directory to return to later
# We use a function so we force return to project directory even if it fails/is killed
run_ABC_seq <- function() { 
  # Admin
  t0 = Sys.time()
  on.exit(setwd(project_dir), add = TRUE)
  set.seed(1)
  
  # Clear any previous ABC step files
  if (dir.exists(output_path)) {
    unlink(file.path(output_path, "*"), recursive = FALSE)
  } else {
    dir.create(output_path, recursive = TRUE)
  }
  setwd(output_path) # Need to set the working directory to the output path for ABC results

  # Run the ABC sequential algorithm
  ABC_seq_res <- ABC_sequential(
      method = abc_params[["method"]],
      #model = run_model, # run single core
      model = run_model_parallel, n_cluster = 12, use_seed = TRUE, # run in parallel
      prior = priors, nb_simul = abc_params[["n_particles"]], 
      summary_stat_target = summ_stat_obs[[1]], # ([[2]] contains the statistics labels)
      alpha = abc_params[["alpha"]], p_acc_min = abc_params[["p_acc_min"]],
      verbose = T)
  
  # Report timings and return the result
  tF = Sys.time()  
  cat("Time taken: ", tF - t0, "\n")
  return(ABC_seq_res)
}
ABC_seq_res <- run_ABC_seq()



# 6. Process the output for easy reading and plotting ------------------------------------------------------
iterations_data <- lapply(ABC_seq_res[["intermediary"]], function(x) {
  colnames(x[["posterior"]])[-1] <- c(names(priors), summ_stat_obs[[2]])
  x[["posterior"]] <- tibble(as.data.frame(x[["posterior"]][,1:5]),
                             stats=asplit((x[["posterior"]][,-(1:5)]),1))
  return(x)
  })
# Remove from ABC_seq to save memory
ABC_seq_res[["intermediary"]] <- NULL

# Add column names to the param and stats matrices
n_particles = nrow(ABC_seq_res[["param"]])
colnames(ABC_seq_res[["param"]]) <- names(priors)
rownames(ABC_seq_res[["param"]]) <- paste0("p",1:n_particles)
colnames(ABC_seq_res[["stats"]]) <- summ_stat_obs[[2]] # Labels for the summary statistic
rownames(ABC_seq_res[["stats"]]) <- paste0("p",1:n_particles)
names(ABC_seq_res[["stats_normalization"]]) <- summ_stat_obs[[2]]
names(ABC_seq_res[["weights"]]) <- paste0("p",1:n_particles)

# Save full version of main results
save(priors, abc_params, obs_data, iterations_data, ABC_seq_res, 
     file=file.path(output_path, "abc_output.RData"))

# Save reduced version of main results for easier load
iterations_data <- lapply(iterations_data, function(x){
                            x[["posterior"]] <- dplyr::select(x[["posterior"]], -stats)
                            return(x)})
save(priors, abc_params, base_params, synthetic_study, obs_data, iterations_data, ABC_seq_res,
     file=file.path(output_path, "abc_reduced_output.RData"))
