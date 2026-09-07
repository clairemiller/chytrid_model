# Libraries and file paths -------------------------------------------------------------
library(Hmisc)
library(tidyverse)
devtools::load_all(".")
source("scripts/figure_formatting.R")

# Load data - experimental
# lab = "exp"
# folder_name = file.path( "seq-abc_exp-results",
#     "easyABC_output-summ_stat_obs_plus_cum_new_inf-alpha0.4-n_particles50000" )

# Load data - synthetic
lab = "syn"
folder_name = file.path("seq-abc_syn-results",
                        "easyABC_output-summ_stat_obs_plus_cum_new_inf-alpha0.4-n_particles50000")


# Output file names
output_filename_log = paste0("figures/",lab,"_distributions_logscale")
output_filename_distributions = paste0("figures/",lab,"_distributions")
output_filename_pairs = paste0("figures/",lab,"_pairs-plot")

# Pre-process -------------------------------------------------------------------
load(file.path("data",folder_name,"abc_reduced_output.RData"))

# Print the number of iterations
cat("Number of iterations: ", length(iterations_data), "\n")

# Convert final iterations into desired tibble format
iterN = tibble::as_tibble(ABC_seq_res[["param"]])
iterN$tab_weight = ABC_seq_res[["weights"]]
# Add beta*alpha 
iterN$alphabeta_sh = iterN$alpha*iterN$beta_sh
iterN$alphabeta_un = iterN$alpha*iterN$beta_un

# Useful to have in long format
iterN_long = pivot_longer(iterN,cols=!tab_weight, names_to = "parameter")

# Calculate summary stats
summ_stats <- iterN_long %>%
  mutate(parameter = factor(parameter, levels = unique(parameter))) |>
  group_by(parameter) |>
  summarise(
    Mean= weighted.mean(x=value, w=tab_weight), 
    Median = wtd.quantile(x=value, w=tab_weight, probs=0.5, normwt = TRUE),
    p2.5 = wtd.quantile(x=value, w=tab_weight, probs=0.025, normwt = TRUE),
    p97.5 = wtd.quantile(x=value, w=tab_weight, probs=0.975, normwt = TRUE))


# ODE system results - for synthetic data study
true_params_syn = lsfitparams[c("beta_sh","beta_un","alpha","omega")]
true_params_syn['alphabeta_sh'] = true_params_syn['alpha']*true_params_syn['beta_sh']
true_params_syn['alphabeta_un'] = true_params_syn['alpha']*true_params_syn['beta_un']
true_params_syn <- data.frame(parameter=names(true_params_syn), 
                          value=true_params_syn, row.names=NULL) %>%
  mutate(parameter = factor(parameter, levels = param_order))



# Table of summary statistics ---------------------------------------------
priors_text <- sapply(priors, function(x) paste0("$U(",x[2],",",x[3],")$"))
summ_stats |>
  mutate(
    Prior = priors_text[parameter],
    parameter = format_labs(parameter, latex=T),
    Mean = round(Mean,digits=3),
    Median = round(Median,digits=3),
    `95%crI`= paste0("(",round(p2.5,digits=3),",",round(p97.5,digits=3),")")) |>
  select(parameter, Prior, Mean, Median, `95%crI`) |>
  mutate(across(everything(), ~ paste0(.x, " &"))) |> # Adds an & for LaTeX table formatting
  knitr::kable(format="simple",
    caption=paste0("Summary statistics for the posteriors (",lab,")"))



# Histogram of posteriors -------------------------------------------------
nbins = 30
prior_bins <- imap_dfr(priors, ~{
  a <- as.numeric(.x[2])
  b <- as.numeric(.x[3])
  tibble(
    parameter = .y,
    xmin = seq(a, b, length.out = nbins + 1)[1:nbins],
    xmax = seq(a, b, length.out = nbins + 1)[2:(nbins + 1)],
    ymin = 0, ymax = (b - a) / nbins
  )
})
  
p_hist <- 
  ggplot(iterN_long) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                colour = "Prior (illustrative)", fill = "Prior (illustrative)"), 
            data = prior_bins, alpha = 0.4, linetype = 2) +
  geom_histogram(aes(x=value, weight = tab_weight,
                     fill="Posterior samples", colour = "Posterior samples"),
                 bins=nbins, boundary = 0, alpha = 0.4) +
  geom_vline(aes(xintercept=Median, linetype = "Median"), data=summ_stats, 
             linewidth = 1) +
  facet_wrap(~parameter, scales="free",
            labeller = as_labeller(format_labs, default = label_parsed)) +
  scale_fill_manual(values = dist_pal) +
  scale_colour_manual(values = dist_pal) +
  scale_linetype_manual(values = 2) +
  labs(y = "Proportion", x="Parameter value", colour = NULL, fill = NULL, linetype=NULL)

# Add true value if it's the synthetic data
if (lab=="syn") {
  p_hist <- p_hist + geom_vline(aes(xintercept=value, linetype="True value"), 
                      linewidth = 1, data=true_params_syn) +
                    scale_linetype_manual(values = c(2,1))
}

# Print to viewer
print(p_hist)



# Histogram - log scale ---------------------------------------------------

p_hist_log <- 
  ggplot(dplyr::filter(iterN_long, parameter %in% names(priors))) +
  geom_histogram(aes(x=value, weight = tab_weight,
                     fill="Posterior samples", colour = "Posterior samples"),
                 bins=nbins, boundary = 0, alpha = 0.4) +
  geom_vline(aes(xintercept=Median, linetype = "Median"), 
             data=dplyr::filter(summ_stats, parameter %in% names(priors)),
             linewidth = 1) +
  facet_wrap(~parameter, scales="free",
             labeller = as_labeller(format_labs, default = label_parsed)) +
  scale_fill_manual(values = dist_pal) +
  scale_colour_manual(values = dist_pal) +
  scale_linetype_manual(values = 2) +
  scale_x_log10() +
  labs(y = "Proportion", x="Parameter value", colour = NULL, fill = NULL, linetype=NULL)

# Add true value if it's the synthetic data
if (lab=="syn") {
  p_hist_log <- p_hist_log + 
    geom_vline(aes(xintercept=value, linetype="true value"), linewidth = 1, 
               data=dplyr::filter(true_params_syn, parameter %in% names(priors))) +
    scale_linetype_manual(values = c(2,1))
}
print(p_hist_log) # Print to viewer




# Pairs plots -------------------------------------------------------------
library(GGally)

# Custom weighted histogram for diagonal
weighted_diag <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(
      aes(weight = weight, y = after_stat(density)),
      bins = 30, fill = "slategrey", colour = "white") +
    custom_theme
}

# Custom contour plot for upper triangle - not weighted
basic_contour <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    scale_fill_gradient(low = "slategray1",high = "slategrey") +
    custom_theme + theme(legend.position = "none")
}

weighted_contour <- function(data, mapping,...) {
  # Calculate the weighted kernel using ks
  xvar <- rlang::as_name(mapping$x)
  yvar <- rlang::as_name(mapping$y)
  df <- data.frame(
    x = data[[xvar]], y = data[[yvar]], w = data[["weight"]])

  w <- df$w * nrow(df) / sum(df$w) # ks::kde wants weights that average to 1
  fit <- ks::kde(x = as.matrix(df[, c("x", "y")]), 
                 w = w, density = T)
  density_df <- expand.grid(x = fit$eval.points[[1]], 
                            y = fit$eval.points[[2]])
  density_df$z <- as.vector(fit$estimate)

  ncontours = 5
  threshold = max(density_df$z, na.rm=T) / 50
  breaks <- c(0, seq(threshold, max(density_df$z, na.rm = TRUE),
                     length.out = ncontours + 1))
  pal <- colorRampPalette(c("slategray1", "slategrey"))(ncontours)
  ggplot(density_df, aes(x, y, z = z)) +
    geom_contour_filled(alpha = 0.9, breaks = breaks) +
    scale_fill_manual( values = c(NA,pal) ) +
    custom_theme + theme(legend.position = "none")
}

# Weighted scatterplot for lower triangle
weighted_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(aes(alpha = weight),
      size = 0.05, colour = "slategray3") +
      scale_alpha_continuous(range = c(0.01, 1)) +
    custom_theme
  }

dat <- cbind(iterN[,names(priors)], weight = iterN$tab_weight)
p <- ggpairs(
  dat,
  columns = 1:length(priors),
  lower = list(continuous = weighted_points),
  upper = list(continuous = weighted_contour),
  diag  = list(continuous = weighted_diag),
  labeller = as_labeller(format_labs, default = label_parsed),
  switch = "both"
) + theme(strip.placement="outside", 
          strip.background = element_blank(),
          panel.spacing = unit(1.5, "lines"))

print(p)



# Effect of vaccination -------------------------------------------------------
vaccination_effect = data.frame(
  eta_P = 1 - iterN$alpha,
  tab_weight = iterN$tab_weight
)
wtd.quantile(x=vaccination_effect$eta_P, w=vaccination_effect$tab_weight, 
             probs=c(0.025,0.5,0.975), normwt = TRUE) |>
  knitr::kable(digits=3, format="simple",
               caption=paste0("Summary statistics for the effect of prior infection (",lab,")"))

library(patchwork)
p_vacc <- ggplot(vaccination_effect) +
  geom_histogram(aes(x=eta_P, weight = tab_weight),
                     fill=dist_pal[2], colour = dist_pal[1],
                 bins=nbins, boundary = 0, alpha = 0.4) +
  geom_vline(aes(xintercept=median(eta_P)), linetype=2) +
  labs(y="Proportion", x=parse(text="eta[P]"))+
  theme(plot.margin = margin(t = 0)) + custom_theme
p_box_vacc <- ggplot(vaccination_effect, aes(x=eta_P, weight = tab_weight)) +
  geom_boxplot(fill = dist_pal[2], colour=dist_pal[1], 
               alpha = 0.4, outlier.shape = NA) + theme_void()

p_box_vacc / p_vacc + plot_layout(heights = c(0.5, 4)) + plot_annotation(title="Effect of prior infection")



# Effect of heating -------------------------------------------------------
heating_effect = data.frame(
  eta_un = 1 - iterN$beta_un/iterN$beta_sh,
  tab_weight = iterN$tab_weight
)
wtd.quantile(x=heating_effect$eta_un, w=heating_effect$tab_weight, 
             probs=c(0.025,0.5,0.975), normwt = TRUE) |>
  knitr::kable(digits=3, format="simple",
               caption=paste0("Summary statistics for the effect of heating (",lab,")"))

xlims <- range(heating_effect$eta_un)
p_heat <- ggplot(heating_effect) +
  geom_histogram(aes(x=eta_un, weight = tab_weight),
                 fill=dist_pal[2], colour = dist_pal[1],
                 bins=nbins, boundary = 0, alpha = 0.4) +
  geom_vline(aes(xintercept=median(eta_un)), linetype=2) +
  labs(y="Proportion", x=parse(text="eta[un]")) +
  xlim(xlims) +
  theme(plot.margin = margin(t = 0)) + custom_theme
p_box_heat <- ggplot(heating_effect, aes(x=eta_un, weight = tab_weight)) +
  geom_boxplot(fill = dist_pal[2], colour=dist_pal[1], 
               alpha = 0.4, outlier.shape = NA) +
  xlim(xlims) + theme_void()

p_box_heat / p_heat + plot_layout(heights = c(0.5, 4)) + 
  plot_annotation(title="Effect of heating")


# Output -----------------------------------------------------------

# Output the histogram
pdf( paste0(output_filename_distributions, ".pdf"), width=11, height=6)
print(p_hist + custom_theme)
dev.off()

# Log scale
pdf(paste0(output_filename_log, ".pdf"), width=11, height=6)
print( p_hist_log + custom_theme )
dev.off()

# Output the pairs plot
png(paste0(output_filename_pairs, ".png"), width=12, height=12, res=150, units="in")
print(p + theme(strip.text = element_text(size = 18)))
dev.off()


