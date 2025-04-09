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

# Plot ordering
param_order <- c("beta_sh", "beta_un", "alpha",
                 "alphabeta_sh", "alphabeta_un", "omega")
compartment_order <- c("S" ,"I1",  "I2",  "I3",  
                       "SV",  "IV1", "IV2", "IV3")

# Fix compartment labels: S/I1-I3 becomes SU/IU1-3
fix_compartment_labels <- function(X) {
  mapping <- c("S"="S[U]",  "I1"="I[U*','*1]",  "I2"="I[U*','*2]",   "I3"="I[U*','*3]",
               "SV"="S[V]", "IV1"="I[V*','*1]", "IV2"= "I[V*','*2]", "IV3"="I[V*','*3]",
               "I"="I[U]", "IV"="I[V]")
  X_new = mapping[X]
  unname(X_new)
}

# Formatting for parameter labels
format_labs <- function(value) {
  lab <- gsub("(beta)_(sh|un)", "\\1[\\2]", value)
  lab <- gsub("(alpha)(beta)", "\\1*\\2", lab)
  return(lab)
}