# Setup ------------------------------------------------------------------
dataset = "exp" # One of "syn" or "exp"

# Check command line arguments if they exist
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  dataset = args[1]
  check1 = assert_that(args %in% c("exp","syn"), 
              msg="Arguments to script must be one of 'exp' or 'syn'.")
}
check2 = assert_that(length(args) <= 1, msg="Too many arguments passed to script.")

# Filepaths
data_filename <- file.path("../data",paste0(dataset,"_data.csv"))
samples_filename <- paste0("../data/",dataset,"_priors_trajectories.RData")
output_filename <- paste0("../data/",dataset,"_summarystats.RData")

# Process data -------------------------------------------------------
# Read in the experimental/synthetic data
expdata <- read.csv(data_filename)
cat("Loaded experimental data:", data_filename, "\n")

# Convert to wide format
expdata <- pivot_wider(expdata, values_fill = 0,
                       names_from="compartment", values_from="N")

# Load the sampled priors and trajectories data
load(samples_filename)
cat("Loaded sampled data:", samples_filename, "\n")

# Function to calculate the summary statistic -------------------------------

# Note: sim_i is one simulation for one experiment (i.e. one environment)
calc_summ_stat <- function(obs, sim_i) {
  # Function for calculating square error for 1 obs and 1 sim
  calc_sq_error <- function(obs_i, sim_i) {
    stopifnot(sim_i$t==obs_i$week)
    ordered_cols <- colnames(sim_i)
    ordered_cols <- ordered_cols[!(ordered_cols %in% c("t","iter_id","shaded"))]
    sq_error = (sim_i[,ordered_cols]-obs_i[,ordered_cols])^2
    return(sq_error)
  }
  # Square error for each mesocosm
  sq_error_arr <- lapply(split(obs, obs$Mesocosm),
                         calc_sq_error, sim_i=sim_i)
  sq_error_arr <- bind_rows(sq_error_arr)
  # Get the sum for each compartment
  sq_error = colSums(sq_error_arr)
  # Return per compartment and full sum
  c(sq_error,"sum_error"=sum(sq_error))
}

# Run summary statistic calculation -----------------------------------------------
# Loop over the iterations in parallel
registerDoParallel(ncores)
tic("ABC statistic calculation time")
error_list <- foreach(i=1:length(simulation_list)) %dopar% {
  sim_i = simulation_list[[i]]
  # Run for each environment for experimental data
  error_shaded = calc_summ_stat(obs=expdata[expdata$shaded,], 
                                sim_i=sim_i[sim_i$shaded,])
  error_unshaded = calc_summ_stat(obs=expdata[!expdata$shaded,], 
                                  sim_i=sim_i[!sim_i$shaded,])
  # Bind together and add the iteration id before returning
  error_i = bind_rows(c(error_shaded,shaded=T), c(error_unshaded,shaded=F))
  error_i$iter_id=sim_i$iter_id[1]
  return(error_i)
}
toc()
per_compartment_error <- bind_rows(error_list)
per_compartment_error$shaded = as.logical(per_compartment_error$shaded)

sample_error = group_by(per_compartment_error, iter_id) %>%
  summarise(sum_error = sum(sum_error))

# Save the data -----------------------------------------------------------
save(priors, sample_error, per_compartment_error, file=output_filename)
cat("Saved summary statistics to", output_filename, "\n")













