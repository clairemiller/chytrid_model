####
# This contains multiple functions to calculate summary statistics for the sequential ABC model.
####

# Choose the summary statistic function 
# (defined in a string here for convenience/labelling, function evaluated at end of this script)
summary_stat_fn_name = "summ_stat_obs_plus_cum_new_inf"
# summary_stat_fn_name = "summ_stat_all_obs"
# summary_stat_fn_name = "summ_stat_prop_I_S" 


# Function to match the simulated data to the observed data ----------------------------------------------------------
match_simulations_to_observations <- function(res_shaded, res_unshaded, obs_data) 
{
  # Get the correct time points for the simulated data
  obs_weeks = unique(obs_data$week)
  stopifnot(obs_weeks == unique(expdata$week))

  sim_shaded = get_timepoints(res_shaded, t_out = obs_weeks)
  sim_unshaded = get_timepoints(res_unshaded, t_out = obs_weeks)
  rownames(sim_shaded) = sim_shaded$t
  rownames(sim_unshaded) = sim_unshaded$t

  # This is just to double check the ordering later on
  sim_shaded$shaded = T
  sim_unshaded$shaded = F

  # Match rows and columns based on shaded column
  ym = rbind(sim_shaded, sim_unshaded)
  n_shaded = nrow(sim_shaded) # Offset the unshaded as it will be the second half of the matrix
  ri <- ifelse(obs_data$shaded, 
               match(obs_data$week, rownames(sim_shaded)),
               match(obs_data$week, rownames(sim_unshaded)) + n_shaded)
  ci <- match(obs_data$compartment, colnames(sim_shaded)) # both should have the same column name order
  stopifnot(sum(colnames(sim_shaded) != colnames(sim_unshaded))==0)
  
  # Confirm the sum_stat_sim and the observed data match
  stopifnot(sum(ym$shaded[ri] != obs_data$shaded)==0)
  stopifnot(sum(ym$t[ri] != obs_data$week)==0)
  stopifnot(sum(colnames(ym)[ci] != obs_data$compartment)==0)
  
  # sim data in correct order for obs data
  y_sim_ordered = ym[cbind(ri,ci)]

  # Return as a data frame
  y_df = data.frame(
    Mesocosm = obs_data$Mesocosm, week=obs_data$week, compartment = obs_data$compartment, N=y_sim_ordered)
  return(y_df)
}


# Summ. Stat. 1. All observations --------------------------------------------------------------------------
summ_stat_all_obs <- function(sim_data_df, obs_data_df)
{
  summ_stat = sim_data_df$N

  # Make labels for the summary statistic
  lab_summ_stat = paste(  paste0("m",sim_data_df$Mesocosm), 
                          paste0("wk",sim_data_df$week), 
                          sim_data_df$compartment,
                          sep=".") 

  # Checks
  stopifnot(length(summ_stat) == nrow(obs_data_df))
  stopifnot(length(summ_stat) == 512)
  
  # Return the summary statistic and labels
  return(list(summ_stat, lab_summ_stat))
}


# Summ. Stat. 2. All observations + cumulative new infections ----------------------------------------------------------------
summ_stat_obs_plus_cum_new_inf <- function(sim_data_df, obs_data_df) {     
  # Calculate the cumulative new infections
  sim_cum_new_inf <- calc_cumulative_new_infections(sim_data_df, unique(obs_data_df$week))
  # Combine with the observed data
  obs_summ_stat_with_labels = summ_stat_all_obs(sim_data_df, obs_data_df)
  summ_stat = c(  obs_summ_stat_with_labels[[1]], 
                  sim_cum_new_inf$cumulative_new_infections)
  # Get the labels for the summary statistic
  lab_cum_new_inf = paste(paste0("m",sim_cum_new_inf$Mesocosm), 
                          paste0("wk",sim_cum_new_inf$week), 
                          "cum_new_inf",
                          sep=".") 
  lab_summ_stat = c(obs_summ_stat_with_labels[[2]], lab_cum_new_inf)

  # Checks
  # stopifnot(length(summ_stat) == 568)

  # Return the summary statistic and labels
  return(list(summ_stat, lab_summ_stat))
}

calc_cumulative_new_infections <- function(df, obs_weeks) {
  # Filter to just be the susceptibles
  df <- df[df$compartment=="S", c("Mesocosm", "week", "N")]
  # Run for each Mesocosm
  cum_inf <- lapply(
    unique(df$Mesocosm), function(meso) {
      df_meso <- df[df$Mesocosm == meso, c("week", "N")]
      # Make sure we have all observation times in the correct order
      obs_weeks = sort(as.numeric(obs_weeks))
      N = setNames(rep(0, length(obs_weeks)), as.character(obs_weeks))
      N[as.character(df_meso$week)] = df_meso$N
      # Now calculate cumulative new infections
      cum_new_infections = cumsum(-diff(N)) # No data point for week=0 due to diff
      # Return as a data frame
      return(data.frame(Mesocosm=meso, 
                        week = as.numeric(names(cum_new_infections)), 
                        cumulative_new_infections = unname(cum_new_infections)))
    }) 
  
  return( do.call(rbind, cum_inf) )
}


# Summ. Stat. 3. Proportion infected and susceptible ----------------------------------------------------------------
summ_stat_prop_I_S <- function(sim_data_df, obs_data_df) {
  
  # Get the total infected and susceptible for each Mesocosm and week
  infected = sim_data_df[grepl("I", sim_data_df$compartment), ]
  infected = dplyr::group_by(infected, Mesocosm, week)
  infected = dplyr::summarise(infected, I=sum(N), .groups = "drop")
  susceptible = sim_data_df[sim_data_df$compartment == "S", c("Mesocosm","week", "N")]
  colnames(susceptible) = gsub("^N$","S",colnames(susceptible))
  total_count = dplyr::group_by(sim_data_df[, c("Mesocosm","week", "N")], Mesocosm, week)
  total_count = dplyr::summarise(total_count, N=sum(N), .groups = "drop")
  combined = merge(merge(infected, susceptible, by=c("Mesocosm","week")), 
                  total_count, by=c("Mesocosm","week"))
  combined = combined[order(combined$Mesocosm, combined$week), ]
  
  # Calculate the proportion infected and susceptible (with no prior infection)
  prop_infected = combined$I / combined$N
  prop_susceptible = combined$S / combined$N
  summ_stat = c(prop_infected, prop_susceptible)
  
  # Make labels for the summary statistic
  label_mesocosm_week <- function(df, main_label) {
    paste(paste0("m",df$Mesocosm),paste0("wk",df$week),main_label, sep=".") }
  lab_infected = label_mesocosm_week(combined, "I")
  lab_susceptible = label_mesocosm_week(combined, "S")
  lab_summ_stat = c(lab_infected, lab_susceptible)
  
  # Some checks
  df_check = expand.grid(Mesocosm = unique(obs_data_df$Mesocosm), week = unique(obs_data_df$week))
  df_check = df_check[order(df_check$Mesocosm, df_check$week), ]
  stopifnot(sum(df_check$Mesocosm != combined$Mesocosm)==0)
  stopifnot(sum(df_check$week != combined$week)==0)
  stopifnot(length(summ_stat) == 128)
  
  # Return the summary statistic and labels
  return(list(summ_stat, lab_summ_stat))
}



# Convert the string name of the summary statistic function to an actual function ------------------------
summary_stat_fn = eval(parse(text = summary_stat_fn_name))
