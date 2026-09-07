library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Colours and theme
env_pal <- c("#2980B9","#C0392B")
dist_pal <- c("#4B165E","#CEA2E1")
custom_theme <- theme_bw() +
  theme(axis.text = element_text(size=unit(16,"pt")),
        axis.title = element_text(size=unit(20,"pt")),
        legend.text = element_text(size=unit(14,"pt")),
        legend.key.height = unit(1,"cm"),
        legend.key.width = unit(1,"cm"),
        legend.key.spacing.y = unit(12,"pt"),
        strip.text = element_text(size=unit(16,"pt")),
        strip.background = element_rect(fill="white"),
        panel.spacing.x = unit(12,"pt"))

# Ordering
param_order <- c("beta_sh", "beta_un", "alpha",
                 "alphabeta_sh", "alphabeta_un", "omega")
compartment_order <- c("S" ,"I1",  "I2",  "I3",  
                       "SV",  "IV1", "IV2", "IV3")

# Fix S/I1-I3 to be SU/IU1-3
fix_compartment_labels <- function(X) {
  mapping <- c("S"="S[N]",  "I1"="I[N*','*1]",  "I2"="I[N*','*2]",   "I3"="I[N*','*3]",
               "SV"="S[P]", "IV1"="I[P*','*1]", "IV2"= "I[P*','*2]", "IV3"="I[P*','*3]",
               "I"="I[N]", "IV"="I[P]")
  X_new = mapping[X]
  unname(X_new)
}
labeller_compartments = as_labeller(fix_compartment_labels, default = label_parsed)

# Format parameter labelling
format_labs <- function(value, latex = F) {
  if (latex) {
    #lab <- gsub("(alpha|beta|omega)", "\\\\\\1", value)
    lab <- paste0("$",value,"$")
  } else {
    lab <- gsub("(beta)_(sh|un)", "\\1[\\2]", value)
    lab <- gsub("(alpha)(beta)", "\\1*\\2", lab)
  }
  return(lab)
}
labeller_parameters = as_labeller(format_labs, default = label_parsed)

