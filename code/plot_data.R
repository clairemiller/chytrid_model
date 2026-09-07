# This script generates the plots of the data used for parameterisation
# - Figure 2 when dataset="exp"
# - Figure S1 when data="syn"

# Setup -------------------------------------------------------------------
dataset = "exp" # One of "syn" or "exp"

# Check command line arguments if they exist
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  dataset = args[1]
  check = assert_that(args %in% c("exp","syn"), 
              msg="Arguments to script must be one of 'exp' or 'syn'.")
}
check = assert_that(length(args) <= 1, msg="Too many arguments passed to script.")

# Filepaths
data_filename <- file.path("../data",paste0(dataset,"_data.csv"))
fig_filename <- file.path("../figures",paste0(dataset,"_data.pdf"))

# Function for plotting ---------------------------------------------------
plot_trajectories <- function(df) {
  figure_df <- mutate(df,
                      experiment = (as.numeric(Mesocosm) %% 4),
                      shaded_lab = ifelse(shaded, "Shaded", "Unshaded"),
                      compartment = factor(compartment,
                                           levels = compartment_order))
  levels(figure_df$compartment) <- fix_compartment_labels(
    levels(figure_df$compartment))
  
  ggplot(figure_df, 
         aes(x=week, y = N, colour=shaded_lab, group=Mesocosm)) +
    geom_point() + geom_line(alpha=0.5) +
    facet_grid(shaded_lab ~ compartment, 
               labeller = label_parsed) +
    guides(colour="none") +
    scale_color_manual(values = env_pal) +
    labs(x="Weeks", y="Num. frogs")
}


# Generate plots ----------------------------------------------------------
# Read in the data
data <- read.csv(data_filename)

# Create the plot
p <- plot_trajectories(data)
# print(p) # Uncomment to print the plots in the gui

# Save to file with custom theme in .Rprofile
pdf(file=fig_filename, width=11, height=6)
print(p + custom_theme)
dev.off()


