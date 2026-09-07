# Libraries and source files
library(dplyr)
library(tidyr)
devtools::load_all(".")

# Experimental setup
N_experiments = 4 # Number of Mesocosms per environment
sample_times = c(0,1,2,4,6,8,10,15) # Observation times

# Set up the parameter sets for shaded and unshaded
parameters[names(lsfitparams)] = lsfitparams
p_shaded = parameters
p_shaded["beta"] = p_shaded["beta_sh"]
p_unshaded = parameters
p_unshaded["beta"] = p_unshaded["beta_un"]

# Set the seed and run the trajectories
set.seed(0)
traj_shaded = run_gillespie(p_shaded, Nsims=N_experiments)
traj_unshaded = run_gillespie(p_unshaded, Nsims=N_experiments)

# Extract the observation time points
traj_shaded = lapply(traj_shaded, get_timepoints, t_out = sample_times)
traj_unshaded = lapply(traj_unshaded, get_timepoints, t_out = sample_times)

# Process into a long dataframe with labelled mesocosms to match expdata.csv
# Columns: Mesocosm,shaded,week,compartment,N
traj_shaded_long <- bind_rows(traj_shaded, .id="Mesocosm") %>%
  mutate(shaded = T) %>% rename("week"="t")
traj_unshaded_long <- bind_rows(traj_unshaded, .id="Mesocosm") %>%
  mutate(Mesocosm = as.numeric(Mesocosm) + N_experiments, # ensure uniqueness
         shaded = F) %>% rename("week"="t")
syndata <- rbind(traj_shaded_long, traj_unshaded_long) %>%
  pivot_longer(cols = !c(Mesocosm, week, shaded),
               names_to="compartment", values_to="N")

# Plot for sanity check
source("scripts/figure_formatting.R")
mutate(syndata, 
  compartment = factor(compartment, levels = compartment_order)) %>%
ggplot( aes(x=week, y = N, colour=shaded, group=Mesocosm) ) +
  geom_point() + geom_line(alpha=0.5) +
  facet_grid(shaded ~ compartment) +
  guides(colour="none") +
  scale_color_manual(values = rev(env_pal)) +
  labs(x="Weeks", y="Num. frogs")

# Save to csv and rda
write.csv(syndata,file="data/syndata_ctmc.csv",quote=F,row.names=F)
usethis::use_data(syndata, overwrite = TRUE)