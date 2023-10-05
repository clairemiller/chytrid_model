# Chytrid Model
Scripts to fit ODE model parameters for the article: Anthony W. Waddle et al., 'Heated homes for hoppers: Hotspot shelters enable frogs to survive chytridiomycosis and stimulate resistance'.
This repository contains 3 files:

1. `process_frog_compartments.R`: a script to process a csv of frog infection timeseries data. This script generates two csv files: *data/expdata_compartments.csv* and *data/expdata_perfrog.csv*,
2. `FitChytridModel.Rmd`: an RMarkdown document which runs the fit on the dataset generated in `process_frog_compartments.R` and generates a html document with the plots and tables presented in the article,
3. `chytrid_model.R`: an R script containing the ODE model, the residuals function, and parameter values, used in the `FitChytridModel.Rmd` document.

These scripts require a csv document (with relative filepath *data/rawdata.csv*) containing the following data columns to run:

- `Mesocosm`: the mesocosm ID
- `shaded`: true/false whether the mesocosm was shaded or unshaded
- `frog_id`: the unique ID of the frog
- `week`: the observation week
- `vaccinated`: true/false whether the frog was vaccinated at week 0
- `infected`: true/false whether the frog was infected at week 0
- `inf_level`: the recorded infection level at this observation
