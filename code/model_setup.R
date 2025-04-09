# Number of prior samples and acceptance rate
N_priors = 10^6
acceptance_p = 0.002
N_accept = round(acceptance_p*N_priors)

# Number of cores to use during ABC analysis
ncores = 4

# Observation times from experimental data (in weeks)
sample_times = c(0,1,2,4,6,8,10,15) 

# Priors definition function
# @param num_samples the number of samples requested
generate_priors <- function(num_samples) {
  data.frame(
    'beta_un' = runif(min=0.0, max=1.0, n=num_samples),
    'beta_sh' = runif(min=0.0, max=1.0, n=num_samples),
    'alpha' = runif(min=0, max=0.5, n=num_samples),
    "omega" = runif(min=0.0, max=10.0, n=num_samples),
    "mu" = runif(min=0.0, max=1.0, n=num_samples)
  )
}
