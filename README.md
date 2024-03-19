# Chytrid Model
Scripts to fit ODE model parameters for the article: Anthony W. Waddle et al., 'Heated homes for hoppers: Hotspot shelters enable frogs to survive chytridiomycosis and stimulate resistance'.
The main document to be run is `FitChytridModel.Rmd`, an RMarkdown document which generates a html document with the figures and tables in the manuscript.
This repository also contains the following files, which are used by `FitChytridModel.Rmd`: 

1. `expdata.csv`: the processed experimental data.
2. `process_frog_compartments.R`: the script that was used to process a raw csv of frog infection timeseries data and generates *expdata.csv*.
3. `chytrid_models_and_residuals`: an R script containing the ODE model, the residuals function, and parameter values.
4. `plotting_functions.R`: an R script with all plotting functions.

To run this function on a different dataset, the raw data csv file to be processed by `process_frog_compartments.R`, needs to contain the following columns:

- `Mesocosm`: the mesocosm ID
- `shaded`: true/false whether the mesocosm was shaded or unshaded
- `frog_id`: a unique frog ID
- `week`: the observation time (in # weeks since start of experiment)
- `vaccinated`: true/false whether the frog was vaccinated at week 0
- `infected`: true/false whether the frog was infected at week 0
- `inf_level`: the recorded infection level at this observation
