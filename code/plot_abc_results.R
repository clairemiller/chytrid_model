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

# Data files
data_filename <- file.path("../data",paste0(dataset,"_data.csv"))
summarystats_filename <- paste0("../data/",dataset,"_summarystats.RData")

# Output filenames
accepted_samples_csv = file.path("../data",paste0(dataset,"_posterior_parameters.csv"))
distributions_filename = paste0("../figures/",dataset,"_distributions.pdf")
log_distributions_filename = paste0("../figures/",dataset,"_log_distributions.pdf")
summary_table_filename = paste0("../figures/",dataset,"_posteriors_statistics.txt")
examples_filename = paste0("../figures/",dataset,"_example_trajectories.pdf")

# Process data for plotting -------------------------------------------------------------------
# Load the expermental data
exp_data <- read.csv(data_filename)

# Load summary statistic data
load(summarystats_filename)
# Add beta*alpha into the priors
priors$alphabeta_sh = priors$alpha*priors$beta_sh
priors$alphabeta_un = priors$alpha*priors$beta_un

# Long version of priors (remove mu)
priors_long = select(priors, -mu) %>%
  pivot_longer(cols = !iter_id) %>%
  mutate(name = factor(name, levels = param_order))

# Process the parameters used for synthetic data generation (remove mu)
syn_params = lsfitparams[names(lsfitparams) != "mu"]
syn_params['alphabeta_sh'] = syn_params['alpha']*syn_params['beta_sh']
syn_params['alphabeta_un'] = syn_params['alpha']*syn_params['beta_un']
syn_params <- data.frame(name=names(syn_params), 
                          value=syn_params, row.names=NULL) %>%
              mutate(name = factor(name, levels = param_order))


# Accepted samples --------------------------------------------------------
# Determine accepted samples
accepted_ids <- sample_error$iter_id[order(sample_error$sum_error)][1:N_accept]
accepted_samples = dplyr::filter(priors, iter_id %in% accepted_ids)

# Save accepted dataframe as csv
write.csv(accepted_samples,file=accepted_samples_csv,
          row.names = F, quote = F)
cat("Saved accepted parameter sets to",accepted_samples_csv, "\n")

# Plot prior and posterior distributions -------------------------------------------------------
make_plot <- function(accepted_df) {
  line_col_mean <- "gray9" 
  line_col_truevalue <- "#3B5998"
  figure_df <- select(accepted_df, -mu) %>%
    pivot_longer(cols = !iter_id) %>%
    mutate(name = factor(name, levels = param_order))
  p <- ggplot(figure_df,
           aes(x=value, y=after_stat(width*density))) + 
      geom_histogram(aes(colour="prior",fill="prior"), bins=30, alpha = 0.4, 
                     data=priors_long) + 
      geom_histogram(aes(colour="posterior",fill="posterior"), alpha = 0.4) +
      facet_wrap(~name, scales="free",
                 labeller = as_labeller(format_labs, default = label_parsed)) +
      scale_fill_manual(values = dist_pal) +
      scale_colour_manual(values = dist_pal) +
      scale_linetype_manual(values = rep("longdash",2)) +
      labs(y = "Proportion of samples", 
           colour = NULL, fill = NULL, linetype = NULL)
  # Add the mean line
  p <- p + stat_summary(fun = mean, geom = "vline", orientation="y",
                     aes(xintercept = after_stat(x), y=0, linetype="mean"),
                     linewidth = 1, colour=line_col_mean)
  # Add true value if it's synthetic
  if (dataset=="syn") {
    p <- p + geom_vline(aes(xintercept=value, linetype="true value"), 
                   linewidth = 1, colour=line_col_truevalue, data=syn_params)
  }
  return(p)
}
p <- make_plot(accepted_samples)

# Save figure
pdf(distributions_filename, width=11, height=6)
p + custom_theme
dev.off()
# Save figure on log scale
pdf(log_distributions_filename, width=11, height=6)
p + scale_x_log10(labels=label_log()) + custom_theme
dev.off()

# Summary statistics table ------------------------------------------------
make_table <- function(accepted_df, caption) {
  select(accepted_df, -mu) %>%
    pivot_longer(cols=!iter_id, names_to="Parameter") %>%
    mutate(Parameter = factor(Parameter, levels = param_order)) %>%
    group_by(Parameter) %>%
    summarise(Mean= mean(value), Median=median(value),
              p2.5 = quantile(value, 0.025),
              p97.5 = quantile(value, 0.975)) %>%
    mutate(Parameter = format_labs(Parameter)) %>%
    knitr::kable(digits=3,caption=caption, format="simple")
}
tab <- make_table(accepted_samples, 
           caption=paste0("Summary statistics for the posterior distributions."))

# Write the table to a text file
cat(tab, file = summary_table_filename, sep="\n")


# Plot posterior examples -------------------------------------------------
# Number of examples to plot
N = 50
N = ifelse(N_accept < N, N_accept, N) # For low sample size testing

# Randomly sample some parameter sets
set.seed(55)
params_to_plot <- accepted_samples[sample(N_accept, size=N),]

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
  return (
    rbind(cbind(traj_shaded, shaded = T), 
          cbind(traj_unshaded, shaded = F)) )
})
trajectories <- bind_rows(trajectories, .id="i")

# Combine into Ss/Is and process for plot
process_data <- function(x) {
  x_long <- mutate(x, I = I1+I2+I3, IV = IV1+IV2+IV3) %>%
    pivot_longer(cols=c("S","SV","I","IV"),
                 names_to="state", values_to="n_frogs") %>%
    mutate(state=factor(state, levels=c("S","I","SV","IV")),
           shaded_lab = ifelse(shaded,"Shaded","Unshaded"))
  levels(x_long$state) <- fix_compartment_labels(levels(x_long$state))
  return(x_long)
}

# Plot the spaghetti plots
traj_long <- process_data(trajectories)
p <- ggplot(traj_long) +
  geom_line(aes(x=t,y=n_frogs,group=i,colour=shaded_lab),alpha=0.2) +
  facet_grid(shaded_lab~state,labeller=label_parsed) + 
  xlim(0,15) +
  labs(x="Weeks", y="Num. frogs", colour=NULL) + 
  scale_colour_manual(values=c(env_pal,"black"),
                      labels=c("Simulated","Simulated","Experimental")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1, size=3) ) )


# Add the experimental data 
exp_data_proc <- pivot_wider(exp_data, 
                             names_from=compartment, values_from=N) %>%
  process_data()
p2 <- p + 
  geom_line(aes(x=week,y=n_frogs,group=Mesocosm, colour="X"), 
            data=exp_data_proc) +
  geom_point(aes(x=week,y=n_frogs,group=Mesocosm, colour="X"), 
             data=exp_data_proc)

# Save to pdf 
pdf(file=examples_filename, width=12, height=6)
print(p2 + custom_theme)
dev.off()
