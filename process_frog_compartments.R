# Processes the raw frog data and converts into compartments
# Outputs from this script are used for parameter fitting in FitChytridMode.Rmd
# S -> I1 -> I2 -> I3 0 -> ...
# SV -> IV1 -> IV2 -> IV3 -> SV ...

# load libraries
library(dplyr)
library(purrr) # For map2
library(ggplot2)
library(tidyr)

# Read in the data, we expect a dataset with the following columns
#   - Mesocosm: the mesocosm ID
#   - shaded: true/false whether the mesocosm was shaded or unshaded
#   - frog_id: the unique ID of the frog
#   - week: the observation week
#   - vaccinated: true/false whether the frog was vaccinated at week 0
#   - infected: true/false whether the frog was infected at week 0
#   - inf_level: the recorded infection level at this observation
rawdata <- read.csv("data/rawdata.csv", check.names = F) %>%
  select(Mesocosm, shaded, week, frog_id,
         vaccinated, infected, inf_level) %>%
  # Remove any frogs with missing infection levels (assumed dead)
  filter(!is.na(inf_level))

# Infection thresholds
thresholds <- list(
  naïve = c( "I1"=10^4,  "I2"=10^5,  "I3"=0 ), 
  vaccinated = c("IV1"=10^3, "IV2"=10^4, "IV3"=0) 
)

# Function to process compartment number
assign_I <- function(inf_level, vaccinated, infected) {
  # Map to the data
  lapply(seq_along(inf_level), function(i) {
    i2.thresh = ifelse(vaccinated[i], 
                       thresholds[["vaccinated"]][["IV2"]],
                       thresholds[["naïve"]][["I2"]])
    
    # Different cases
    if (inf_level[i] == 0) {
      compartment = 0
    } else if (inf_level[i] > i2.thresh) {
      compartment = 2
    } else if (i==1) {
      compartment = 1 # All not infected should be picked up by first condition
    }
    else {
      compartment = c(1,3)[which.min(c(inf_level[i-1], inf_level[i]))]
    }
    return(compartment)
  }) %>% unlist()
}

# Calculate the compartment sequence for each frog
output_perfrog <- rawdata %>% 
  arrange(frog_id,week) %>%
  # First compartment assignment, based on infection level and vaccinated or not
  group_by(frog_id) %>%
  mutate(compartment_i = assign_I(inf_level,vaccinated,infected)) %>%
  # Check for previous infection
  group_by(frog_id) %>%
  mutate(prev_infection = cumsum(c(0,diff(compartment_i))<0)>0) %>%
  # Re-determine infection level using previous infection and vaccinated
  mutate(compartment_i = assign_I(inf_level,vaccinated | prev_infection, infected)) %>%
  # Create label from compartment and vaccinated/previous infection
  mutate(compartment = paste0(ifelse(compartment_i==0,"S","I"),
                              ifelse(vaccinated | prev_infection,"V",""),
                              ifelse(compartment_i==0,"",compartment_i)))

# Collate the per frog data to compartmental data for each mesocosm
output_permesocosm <- output_perfrog %>%
  group_by(Mesocosm,shaded,week,compartment) %>%
  summarise(N=n(),.groups="drop") %>%
  arrange(Mesocosm,week)

# Output both datasets to CSV
write.csv(x=output_permesocosm,file="data/expdata_compartments.csv",
          quote=F, row.names=F)

output_perfrog %>%
  arrange(Mesocosm,frog_id,week) %>%
  relocate(frog_id,.after=shaded) %>%
  write.csv(x=.,file="data/expdata_perfrog.csv",quote=F,row.names=F)

