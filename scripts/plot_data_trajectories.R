# Load functions
devtools::load_all(".")
source("scripts/figure_formatting.R")

# Function for plotting ---------------------------------------------------
plot_trajectories <- function(df) {
  figure_df <- mutate(df,
                      experiment = (as.numeric(Mesocosm) %% 4),
                      shaded_lab = ifelse(shaded, "Shaded", "Unshaded"),
                      compartment = factor(compartment,
                                           levels = compartment_order))
  
  ggplot(figure_df, aes(x=week, y = N, colour=shaded_lab, group=Mesocosm)) +
    geom_point() + geom_line(alpha=0.5) +
    facet_grid(shaded_lab ~ compartment, 
               labeller = labeller(compartment = labeller_compartments)) +
    guides(colour="none") +
    scale_color_manual(values = env_pal) +
    labs(x="Weeks", y="Num. frogs")
}

# Experimental ------------------------------------------------------------
exp_fig_filename <- "figures/exp_data.pdf"

# Plot the trajectories
exp_p <- plot_trajectories(expdata)
print(exp_p)

# Save to file
pdf(file=exp_fig_filename, width=11, height=6)
print(exp_p + custom_theme)
dev.off()


# Synthetic ---------------------------------------------------------------
syn_fig_filename <- "figures/syn_data.pdf"

# Plot the trajectories
syn_p <- plot_trajectories(syndata)
print(syn_p)

# Plot the trajectories
pdf(file=syn_fig_filename, width=11, height=6)
print(syn_p + custom_theme)
dev.off()

