# Paper

This repository contains the data and code required to reproduce the results for the paper: 'Sunlight-heated refugia protect frogs from chytridiomycosis: a mathematical modelling study', https://doi.org/10.48550/arXiv.2503.06846. 

## Abstract

The fungal disease Chytridiomycosis poses a threat to frog populations worldwide. It has driven over 90 amphibian species to extinction and severely affected hundreds more. Difficulties in disease management have shown a need for novel conservation approaches.

We present a novel mathematical model for chytridiomycosis transmission in frogs that includes the natural history of infection, to test the hypothesis that sunlight-heated refugia reduce transmission. This model was fit using approximate Bayesian computation to experimental data where frogs were grouped into sunlight-heated or shaded refugia cohorts.

Results show a 40% reduction in infection due to sunlight-heating of refugia. Frogs that were infected and recovered had a reduction in susceptibility of approximately 97% compared to naive frogs. Our model offers insight into using sunlight-heated refugia to reduce chytridiomycosis prevalence. Importantly, it is the first step in determining the necessary level of refugia in the landscape for frog population recovery and population sustainability.

# Chytrid Model

## Repository Structure

The structure and files in this repository are as follows:

``` bash
├── Dockerfile
├── README.md
├── run_experimental_abc.sh
├── run_synthetic_abc.sh
├── code
│   ├── .Rprofile
│   ├── calc_summary_stats.R
│   ├── chytrid_model.Rproj
│   ├── fns_stochastic_model.R
│   ├── generate_synthetic_data.R
│   ├── model_setup.R
│   ├── plot_abc_results.R
│   ├── plot_data.R
│   ├── plot_formatting.R
│   └── simulate_samples.R
├── data
│   ├── exp_data.csv
│   └── syn_data.csv
└── figures
```

## Running the Model and Reproducing Paper Results

To reproduce the results and figures from the paper, follow these steps:

1.  **Data Storage**: All intermediate data generated at each step is stored in the `data` directory, while the resulting figures and tables are saved in the `figures` directory.

2.  **Modifying Parameters**: To adjust the number of samples generated or the acceptance rate, update the relevant values in the `model_setup.R` file.

3.  **Running the Workflow**:

    i.  The complete workflow for processing the experimental and synthetic data can be executed using the bash scripts in the home directory: `bash run_experimental_data.sh` or `bash run_synthetic_data.sh` respectively.

    ii. Otherwise, following the instructions below for a step-by-step analysis. These use the notation `exp` for experimental data (main manuscript results) and `syn` for synthetic data (Supplementary S2). Also note the use of `.Rprofile` for package and loading of source R file parameters and functions (including `model_setup.R`).

4.  **Command Execution**: All commands are executed from within the `code` directory.

5.  **Estimated run times**: $10^6$ samples on only 4 cores is expected to take around 3-4 hours and we recommend using HPC for the this sample count. We recommend testing the code using $10^4$ samples, which only takes approximately 5 minutes on 4 cores.

### Main Manuscript (Experimental Data)

1.  **Plot Experimental Data:**

    ``` bash
    Rscript plot_data.R 'exp'
    ```

2.  **Run ABC Simulation:**

    i.  Generate trajectories for all prior samples (saves data to `exp_priors_trajectories.RData`):

        ``` bash
        Rscript simulate_samples.R 'exp'
        ```

    ii. Calculate summary statistics for all samples (saves data to `exp_summarystats.RData`):

        ``` bash
        Rscript calc_summary_stats.R 'exp'
        ```

3.  **Plot prior and posterior distributions (Figure 3(a)), posterior examples (Figure 3(b)), and extract Summary Statistics (Table 3) :**

    ``` bash
    Rscript plot_abc_results.R 'exp'
    ```

### Synthetic Data (Supplementary S2)

1.  **Regenerate Synthetic Data:** (not required, synthetic data used in manuscript is included in `data/syn_data.csv`.)

    ``` bash
    Rscript generate_synthetic_data.R
    ```

2.  **Run Pipeline Described for Main Manuscript, but Replace 'exp' with 'syn':**

    ``` bash
    Rscript plot_data.R 'syn'
    Rscript simulate_samples.R 'syn'
    Rscript calc_summary_stats.R 'syn'
    Rscript plot_abc_results.R 'syn'
    ```

## R Package Requirements

The following R libraries are required to run the code (specific versions used given in parenthesis):

-   `ggplot2` (version 3.5.1)
-   `tidyr` (version 1.3.0)
-   `dplyr` (version 1.1.1)
-   `scales` (version 1.3.0)
-   `GillespieSSA` (version 0.6.2)
-   `assertthat` (version 0.2.1)
-   `doParallel` (version 1.0.17)
-   `tictoc` (version 1.2.1)
-   `knitr` (version 1.42)

To install any missing packages, you can use the following R code:

``` r
packages <- c("ggplot2", "tidyr", "dplyr", "scales", "GillespieSSA", "assertthat", "doParallel", "tictoc", "knitr")

installed_packages <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
}
```

## Running Through Docker

We have also included a Dockerfile to run this code using a Docker container without the requirement for R and library installation. For instructions on installing Docker, see: <https://www.docker.com/get-started/>.

1.  **Build the Docker Image:** The docker image must be built before containers can be run (only needs to be done once). This is done from the root directory of the project (the folder containing the Dockerfile).

    ``` bash
    docker build -t chytrid_model .
    ```

Scripts can be run individually through the containers as detailed above. Alternatively, the easiest way to generate the data and plots for each study is as follows:

2.  **Start the Docker Container:**

    ``` bash
    docker run --volume ./:/project -it --rm chytrid_model
    ```

3.  **Run the Experimental Data Study:**

    ``` bash
    bash run_experimental_data.sh
    ```

4.  **Run the Synthetic Data Study:**

    ``` bash
    bash run_synthetic_data.sh
    ```

5.  **Exit and Stop the Container When Finished:**

    ``` bash
    exit
    ```
