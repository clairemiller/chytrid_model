
# Set the seed and synthetic data output
set.seed(0)
N_experiments = 4 # Number of Mesocosms per environment
output_filename <- "../data/syn_data.csv"

# Set up the parameter sets for shaded and unshaded environments
parameters[names(lsfitparams)] = lsfitparams
p_shaded = parameters
p_shaded["beta"] = p_shaded["beta_sh"]
p_unshaded = parameters
p_unshaded["beta"] = p_unshaded["beta_un"]

# Run the trajectories
traj_shaded = run_gillespie(p_shaded, Nsims=N_experiments)
traj_unshaded = run_gillespie(p_unshaded, Nsims=N_experiments)

# Extract the time points
traj_shaded = lapply(traj_shaded, get_timepoints, t_out = sample_times)
traj_unshaded = lapply(traj_unshaded, get_timepoints, t_out = sample_times)

# Process into a long dataframe with labelled mesocosms to match expdata.csv
# Columns: Mesocosm,shaded,week,compartment,N
traj_shaded_long <- bind_rows(traj_shaded, .id="Mesocosm") %>%
  mutate(shaded = T) %>% rename("week"="t")
traj_unshaded_long <- bind_rows(traj_unshaded, .id="Mesocosm") %>%
  mutate(Mesocosm = as.numeric(Mesocosm) + N_experiments, # ensure uniqueness
         shaded = F) %>% rename("week"="t")
output_df <- rbind(traj_shaded_long, traj_unshaded_long) %>%
  pivot_longer(cols = !c(Mesocosm, week, shaded),
               names_to="compartment", values_to="N")

# Save to file
write.csv(output_df,file=output_filename,quote=F,row.names=F)
