source("chytrid_models_and_residuals.R")
expdata <- expdata <- read.csv("data/expdata.csv", check.names = F) %>%
    mutate(Mesocosm = as.factor(Mesocosm)) %>%
    pivot_wider(names_from = compartment, values_from = N) %>%
    mutate(across(-c(Mesocosm, shaded, week), ~ ifelse(is.na(.x), 0, .x))) %>%
    rowwise() %>%
    mutate(N = sum(c_across(-c(Mesocosm, shaded, week)))) %>%
    ungroup()

obs.long <- read.csv("data/expdata.csv")
# Change the time column from week to t to match function nomenclature
obs <- rename(obs.long, "t" = "week") %>%
    # We want to match the wide format with one time column 't' and one column per compartment
    pivot_wider(names_from = "compartment", values_from = "N", values_fill = 0) %>%
    select(-Mesocosm)
# You'll need to source some stuff from RunABC-HPC.R too

omega_sweep <- seq(0, 5, by=0.001)

test_parameters <- readRDS("data/test_parameters.rds")

omega_residuals <- lapply(cli::cli_progress_along(omega_sweep), function(o) {
    
    test_parameters[["omega"]] <- omega_sweep[o]

    tibble(
        omega = omega_sweep[o],
        residual = sum(calc_residuals(test_parameters, expdata)^2)
    )
}) %>% bind_rows() %>% mutate(type = "ode")

omega_gillespie_residuals <- lapply(cli::cli_progress_along(omega_sweep), function(o) {
    test_parameters[["omega"]] <- omega_sweep[o]

    gill <- runModelBothEnv(params.fit = test_parameters |> unlist(), Tvec = unique(obs$t))

    
    tibble(
        omega = omega_sweep[o],
        residual = calcSummaryStatisticsBothEnv(obs, gill)
    )
}) %>% bind_rows() %>% mutate(type = "gillepsie")

rbind(omega_residuals, omega_gillespie_residuals %>% mutate(residual = 2*residual)) %>%
ggplot(aes(x = omega, y = residual, colour = type)) +
    geom_line()
