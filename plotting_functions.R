# Plotting functions for chytrid model RMarkdown document FitChytridModel.Rmd
# This file includes:
# 1. Function to plot different predictions against experimental data
# 2. Function to run stochastic simulations and plot against deterministic prediction
# Note: this script relies on functions and variables in chytrid_models_and_residuals.R

# Claire Miller & Jen Flegg  & Trish Campbell Feb 2024

# 1. Plot predictions against experimental data ------------------------------
# Plots the timeseries of the compartments for the given 
# ode parameters against the experimental data
# inputs: 
#   - parameters_list (list of named numeric vectors)
#   - compdata (data frame of experimental data)
#   - line_lab_arr (character vector with labels for list items in parameters_list)
plot_compartments <- function(parameters_list, compdata, 
                              line_lab_arr = c("Prediction")) {
  # Function to convert labels to math notation
  lab_fn <- function(labels) {
    gsub("(\\a)|(V*\\d)|(V)","\\1[\\2\\3]",labels)
  }
  # Prediction given parameters
  prediction_list <- setNames(parameters_list, line_lab_arr) %>%
    lapply(function(parameters) {
      ode(
        func=chytrid.ode.model,
        y=x0,
        times=seq(0,max(compdata$week),length.out=100),
        parms=parameters
      ) %>% as_tibble() %>%
        mutate(across(everything(),as.numeric)) %>%
        pivot_longer(cols=-time) %>%
        filter(name != "N") %>%
        mutate(name = lab_fn(name))
    })
  prediction = bind_rows(prediction_list,.id="label")
  # Process experimental data into long format
  compdata_long <- select(compdata,-N) %>%
    pivot_longer(cols=!c(Mesocosm,shaded,week)) %>%
    filter(name != "N") %>%
    rename(time="week") %>%
    mutate(name = lab_fn(name))
  # Plot the dats
  p <- ggplot(compdata_long,aes(x=time,y=value)) + 
    geom_line(aes(group=Mesocosm),alpha=0.5) +
    geom_point(aes(group=Mesocosm,colour="Data"),alpha=0.7) +
    facet_wrap(~name, labeller="label_parsed") +
    labs(x="Time (weeks)",y="Number of frogs",colour=NULL)+
    theme(legend.position=c(0.85,0.15))
  # Add the predictions and set colours/linetypes
  p <- p + 
    geom_line(aes(colour=label, linetype = label),size=1,data=prediction) +
    scale_colour_manual(values=c("black","#1B9E77","darkorange")) +
    scale_linetype_manual(values = c(1,5), guide = "none")
  
  return(p)
}

# 2. Plot Gillespie trajectories for given params ------------------------------
# Plots the trajectories of the compartments for the given parameters against
# the deterministic prediction.
# inputs: 
#   - parms (named numeric vectors)
plot_gillespie <- function(parms) {
  # Run the trajectories using function in chytrid_models_and_residuals.R
  trajectories_list <- run_gillespie(parms)
  # Process into 1 data frame for plotting
  trajectories <- lapply(trajectories_list, function(x) {
    x$data %>%
      data.frame() %>%
      pivot_longer(!t, names_to = "compartment", values_to = "frogcount")}) %>%
    bind_rows(.id="run") %>%
    mutate(compartment = gsub("(^[A-Z])(.*)","\\1[\\2]",compartment))
  # Plot trajectories
  p <- ggplot(trajectories, aes(x=t,y=frogcount)) +
    geom_line(aes(group=run, colour="Stochastic"), alpha=0.4) + 
    facet_wrap(~compartment, labeller=label_parsed) +
    labs(color=NULL, x="Time (weeks)",y="Number of frogs")
  
  # Run the deterministic prediction and process for plotting
  prediction <- ode(
    func=chytrid.ode.model,
    y=x0,
    times=seq(0,max(expdata$week),length.out=100),
    parms=parms
  ) %>% as_tibble() %>%
    mutate(across(everything(),as.numeric)) %>%
    select(-N) %>%
    pivot_longer(cols=-time, names_to="compartment", values_to="frogcount") %>%
    rename(t=time)  %>%
    mutate(compartment = gsub("(^[A-Z])(.*)","\\1[\\2]",compartment))
  # Add to trajectory plot
  p <- p + geom_line(aes(colour="Deterministic"), data=prediction, 
                     size=1) +
    scale_colour_manual(values=c("#1B9E77", "black")) +
    theme(legend.position=c(0.85,0.15))
  
  # Format plot and return
  p <- p + scale_y_continuous(
    breaks = seq(0, 16, by = 4),
    limits = c(0,16),
    expand = expand_scale(mult = c(0,  0.05))
  ) + coord_cartesian(xlim = c(0,  15), expand=0)
  return(p)
}