# Libraries and data load ---------------------------------------------------------------
devtools::load_all(".")
source("scripts/figure_formatting.R")

# Load data - experimental
lab = "exp"
folder_name = file.path( "seq-abc_exp-results",
    "easyABC_output-summ_stat_obs_plus_cum_new_inf-alpha0.4-n_particles50000" )

# Load data - synthetic
# lab = "syn"
# folder_name = file.path("seq-abc-syn-results",
#                         "easyABC_output-summ_stat_obs_plus_cum_new_inf-alpha0.4-n_particles50000")

# Load the data
load(file.path("data",folder_name,"abc_reduced_output.RData"))

# Output file names
fig_filename = paste0("figures/",lab,"_trajectories")



# Select samples and simulate ----------------------------------------------------------
# Number of examples to plot
N = 50
set.seed(55)

# Process into data frame
posterior_df = tibble::as_tibble(ABC_seq_res[["param"]])
posterior_df$weight = ABC_seq_res[["weights"]]


# Randomly sample some parameter sets (using weights)
random_particle_ids = sample(nrow(posterior_df), size=N, 
                             prob=posterior_df$weight)
params_to_plot <- posterior_df[random_particle_ids,]

# Run Gillespie
trajectories <- lapply(1:N, function(i) {
  p_new <- unlist(params_to_plot[i,])
  parameters[names(p_new)] = p_new
  p_shaded = parameters
  p_shaded[["beta"]] = p_shaded[["beta_sh"]]
  p_unshaded = parameters
  p_unshaded[["beta"]] = p_unshaded[["beta_un"]]
  traj_shaded = run_gillespie(p_shaded, Nsims=1)[[1]][["data"]] %>%
    as.data.frame()
  traj_unshaded = run_gillespie(p_unshaded, Nsims=1)[[1]][["data"]] %>%
    as.data.frame()
  # Check for appending final time at end of data frame
  traj_shaded = distinct(traj_shaded)
  traj_unshaded = distinct(traj_unshaded)
  return (
    rbind(cbind(traj_shaded, shaded = T), 
        cbind(traj_unshaded, shaded = F)) )
}) |> 
  bind_rows(.id="Mesocosm")  |>
  rename("week"="t") 



# Process data for plotting -----------------------------------------------

# Combine into Ss/Is
process_data <- function(x) {
  x |> 
    mutate(vaccinated = ifelse(grepl("V",compartment),"[P]","[N]"),
           compartmentSI = gsub("(I|S)(.*)","\\1",compartment),
           shaded = ifelse(shaded, "Shaded", "Unshaded")) |> 
    mutate(compartmentSINV = paste0(compartmentSI,vaccinated)) |>
    group_by(Mesocosm, shaded, week, compartmentSINV) |>
    summarise(N = sum(N), .groups = "drop") |>
    mutate(compartmentSINV = factor(compartmentSINV, 
                                    levels = unique(compartmentSINV)[c(3,1,4,2)]))
}
# Trajectories
traj_proc <- pivot_longer(trajectories, cols=!c(Mesocosm,week,shaded), 
                          names_to = "compartment", values_to="N") |>
              process_data()
# Experimental data
weeks = unique(obs_data$week)
obs_data_proc <- process_data(obs_data) |>
# Fill in missing weeks with 0s
    group_by(Mesocosm, compartmentSINV, shaded) |>
    complete(week = weeks, fill = list(N = 0)) 



# Create the spaghetti plots -----------------------------------------------
p <- ggplot(traj_proc) +
  geom_line(aes(x=week,y=N,group=Mesocosm,colour=shaded),alpha=0.2) +
  facet_grid(shaded~compartmentSINV, labeller = labeller(.cols = label_parsed)) + 
  coord_cartesian(xlim = c(0, 15),expand = 0) +
  labs(x="Weeks", y="Num. frogs", colour=NULL) + 
  scale_colour_manual(values=c(env_pal,"black"),
                      labels=c("Simulated","Simulated",
                               ifelse(lab=="exp","Experimental","Synthetic"))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1, size=3) ) )
# Add experimental data
p <- p + 
  geom_line(aes(x=week,y=N,group=Mesocosm, colour="X"), 
            data=obs_data_proc) +
  geom_point(aes(x=week,y=N,group=Mesocosm, colour="X"), 
             data=obs_data_proc)
print(p)



# Save to pdf -------------------------------------------------------------
pdf(file=paste0(fig_filename,".pdf"), width=11, height=6)
  print(p + custom_theme)
dev.off()

