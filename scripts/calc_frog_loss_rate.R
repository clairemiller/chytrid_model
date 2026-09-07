
library(dplyr)
devtools::load_all(".")

# Function for the loss rate
calc_frog_loss_rate <- function(df) {
    Tf = max(df$week) # Final observation time
    stopifnot(Tf == 15)
    F_df = filter( df, week == 0 | week == Tf ) |> 
            group_by(Mesocosm, week) |> summarise(F=sum(N), .groups="drop")
    F0 = mean( F_df$F[F_df$week == 0] ) # Average initial frog count per mesocosm
    stopifnot(F0 == 20)
    Fbar_Tf = mean( F_df$F[F_df$week == Tf] ) # Average final frog count per mesocosm
    mu = -log( Fbar_Tf / F0 ) / Tf # Loss rate
    return( mu )
}

# Experimental data
mu_exp = calc_frog_loss_rate(expdata)
cat("Loss rate for experimental data: ", mu_exp, "\n")

# Synthetic data
mu_syn = calc_frog_loss_rate(syndata)
cat("Loss rate for synthetic data: ", mu_syn, "\n")
