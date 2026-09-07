cd code

printf '%*s\n' 80 '' | tr ' ' '-'
echo Plotting experimental data....
Rscript plot_data.R 'exp'

printf '%*s\n' 80 '' | tr ' ' '-'
echo Running abc...
Rscript simulate_samples.R 'exp'
Rscript calc_summary_stats.R 'exp'

printf '%*s\n' 80 '' | tr ' ' '-'
echo Plotting results....
Rscript plot_abc_results.R 'exp'

cd ..