# Paper

This repository contains the data and code required to reproduce the results for the paper: 'Estimating the impact of prior infection and sunlight-heated refugia on chytridiomycosis in frogs: a mathematical modelling study'. 

## Abstract

The fungal disease chytridiomycosis has driven over 90 amphibian species to extinction, and severely affected hundreds more. Difficulties in disease management have shown a need for novel conservation approaches.
We present a novel mathematical model for chytridiomycosis transmission in frogs, that includes the natural history of infection, to test the hypothesis that sunlight-heated refugia reduce transmission. The model is fit using approximate Bayesian computation to previous experimental data where a cohort of frogs with varying pre-existing exposure to infection were provided access to either sunlight-heated or shaded refugia. Using our model, we estimate the extent to which prior chytridiomycosis infection protects against subsequent infection and quantify the effect of sunlight-heating of refugia.
Results estimate a 46\% (95\% CrI (-125\%, 87\%)) reduction in infection due to sunlight-heating of refugia.
While our estimates suggest substantial impacts of sunlight-heated refugia on disease transmission, more data is required to resolve uncertainties. Frogs that were infected and recovered had an estimated reduction in susceptibility of 73\% (95\% CrI (27\%, 94\%)) compared to frogs with no prior infection. 
This study contributes to the evidence base for the use of sunlight-heated refugia as part of conservation strategies, and demonstrates the use of mathematical modelling to inform the implementation of habitat-based interventions for amphibian population recovery and sustainability.

# Running the model

## Repository Structure

The structure and files in this repository are as follows:
- The `R` directory contains the core model definition and simulation code and ABC summary statistic methods.
- The `scripts` directory contains all the scripts to run the estimation and plot the results. 
- The `data` directory contains the processed experimental, and is where all estimation results are stored.
- The `figures` directory is where any generated figures are saved.

To switch between synthetic and experimental:
- In `run_seq_abc.R` set the boolean flag `synthetic_study` at the beginning of the file 
- In `plot_posterior_distributions.R`and `plot_posterior_examples.R` add 'syn' or 'exp' as a command line argument when running the script, or adjust the appropriate commented default for `lab` at the beginning of the relevant file.

The directory tree below details the files included in this repository and annotates their use. 

``` bash
├── README.md
├── chytrid-transmission-model.Rproj # R project file
├── DESCRIPTION # Description folder for devtools project setup
├── scripts
│   ├── plot_posterior_examples.R # Plot the example trajectories given the posterior particle set
│   ├── plot_posterior_distributions.R # Plot the posterior distributions
│   ├── run_seq_abc.R # Run the estimation
│   ├── generate_synthetic_data_ctmc.R # Generate the synthetic data for running the simulation estimation study
│   ├── plot_data_trajectories.R # Plot the trajectories of the experimental and synthetic data (generated in generate_synthetic_data_ctmc)
│   └── calc_frog_loss_rate.R # Calculate the loss rate of the frogs for each setup
├── R
│   ├── ode_model.R # The ODE model definition (to determine parameters for the synthetic study)
│   ├── stochastic_model.R # Model and parameter definition for the CTMC model
│   ├── seq_abc_summ_stats.R # Definition of summary statistic for the ABC estimator
│   ├── gillespie_run_fns.R # Functions to run the CTMC model using the gillespie algorithm
│   └── figure_formatting.R # Plot and label formatting
├── data
│   ├── expdata.rda # Project data - experimental data from Waddle et al.
│   ├── syndata.rda # Project data - synthetic data generated using `generate_synthetic_data_ctmc.R`
│   ├── expdata_filled.csv # CSV version of expdata.rda
│   └── syndata_ctmc.csv # CSV version of syndata.rda
└── figures
```


## Running the Model and Reproducing Paper Results

Considerations for running the paper results:

1.  **Data Storage**: All intermediate data generated at each step is stored in the `data` directory, while the resulting figures are saved in the `figures` directory. All intermediate results are stored and require on the order of **15—20 GB of storage** per estimation. This can be reduced by changing the `verbose = T` flag in `run_seq_abc.R` (line 118), but note this may break the data saving and subsequent plotting code.

2.  **Estimated run times**: The setup described in the paper, on 12 cores (Mac Studio M4 Max), takes between **4 and 7 days**.

3.  **Modifying Parameters**: To adjust priors or algorithm hyperparameters, update the relevant values in the `run_seq_abc.R` file. To adjust other model parameters, update the relevant values in `stochastic_model.R`.

4.  **Script Execution**: All scripts are setup to be executed from within the root project directory.

### Model run workflow
*The same workflow is used for both the simulation estimation study and the main estimation results.*

0. **Regenerate synthetic estimation study data (if desired)**
    ```bash
    Rscript generate_synthetic_data_ctmc.R
    ```

1.  **Plot observations:**
    
    This will plot both synthetic and experimental data. Figures are saved to figures directory with prefix `syn_` or `exp_` for simulation estimation and main estimation respectively.

    ```bash
    Rscript plot_data_trajectories.R
    ```

2.  **Run ABC Simulation:**

    i. Set boolean `synthetic_study`in `run_seq_abc.R` (`TRUE` for simulation estimation, `FALSE` for main estimation). 

    ii. By default, data is stored under a subdirectory labelled by algorithm hyperparameters within directories `data/seq-abc_syn-results` and `data/seq-abc_exp-results` for the simulation estimation and main estimation respectively.
    
    ```bash
    Rscript run_seq_abc.R
    ```

3.  **Plot prior and posterior distributions, posterior examples, and extract Summary Statistics:**

    i. Use command line argument 'syn' or 'exp' to switch between the simulation estimation or main estimation data.
    
    ii. Figures are saved to figures directory with prefix `syn_` or `exp_` for simulation estimation and main estimation respectively.

    ``` bash
    Rscript plot_posterior_distributions.R 'syn'
    Rscript plot_posterior_examples.R 'syn'
    ```



## R Package Requirements

The following R libraries are required to run the code (specific versions used given in parenthesis):

-   `tidyverse` (version 2.0.0)
-   `scales` (version 1.4.0)
-   `GillespieSSA` (version 0.6.2)
-   `devtools`(version 2.5.2)
-   `EasyABC` (version 1.6)
-   `Hmisc` (version 5.2.6)
-   `GGally`(version 2.4.0)

To install any missing packages, you can use the following R code:

``` r
packages <- c("tidyverse", "scales", "GillespieSSA", "devtools", "EasyABC", "Hmisc", "GGally")

installed_packages <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
}
```


