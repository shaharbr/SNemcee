# SNemcee

**Dependencies:** <br>
Python libraries: <br>
Corner <br>
Emcee <br>
Pandas <br>
Numpy <br>
Datetime <br>
scipy <br>


**Input data:** <br>
Supernova data files (observations) should be stored in the directory ‘/results’, with the following names: <br>
<SN_name>_veloc <br>
<SN_name>_mag <br>
<SN_name>_lum <br>
(See example files for expected format) <br>

In addition, there should be a table in ‘/results’ that contains the distance (in Mpc) <br>

**Running SNemcee:** <br>
To run the SNemcee program on your SN of choice, use command: <br>

python3 run_mcmc_snec_for_all_SNe.py -S <SN name> -s <number of steps> -w <number of walkers> -o <output directory> <br>

By default, it will run with the following variations of mcmc_snec (can be changed by modifying the code): <br>
<SN_name>_lum_csmFalse_normalizedFalse_TreshLumFalse <br>
<SN_name>_lum_csmTrue_normalizedFalse_TreshLumFalse <br>
<SN_name>_lum-veloc_csmTrue_normalizedFalse_TreshLumFalse_TreshVelocTrue <br>
<SN_name>_lum-veloc_csmTrue_normalizedFalse_TreshLumTrue_TreshVelocTrue <br>
<SN_name>_lum-veloc_csmTrue_normalizedTrue_TreshLumFalse_TreshVelocTrue <br>
<SN_name>_lum-veloc_csmTrue_normalizedTrue_TreshLumTrue_TreshVelocTrue <br>
<SN_name>_mag-veloc_csmTrue_normalizedFalse_TreshMagTrue_TreshVelocTrue <br>
<SN_name>_mag-veloc_csmTrue_normalizedTrue_TreshMagTrue_TreshVelocTrue <br>

**SNemcee outputs:** <br>
Output files will be stored in ‘/mcmc/mcmc_results/<output directory>/<number of walkers>’. <br>
Output files will include summary figures: <br>
Corner plot of posterior distributions for each parameter (all walkers, and steps after burn-in) <br>
Fits of final models to the input observation data for that SN (all walkers, and steps after burn-in) <br>
Chain plots for each parameter <br>
Flat_sampler.csv: <br>
MCMC sampler output table - a 2D table, where columns are the parameters, and rows are walkers and steps (by the hierarchy of: steps, and within them walkers) <br>
final_results.csv: <br>
Summary statistics for the final posterior distributions for each parameter: mean, 16% percentile and 85% percentile. <br>
run_parameters.csv: <br>
Record of the input arguments used in that run <br>

**Running script for comparison plots:** <br>
After running run_mcmc_snec_for_all_SNe.py, you can use plot_snec_fits_wrapper.py to plot paper-ready figures - for example, to compare the effect of different likelihood functions (different fitting schemes), or the inclusion of CSM, on the posterior distributions: <br>

python3 plot_snec_fits_wrapper.py -S <SN name> -s <number of steps> -t <type of figure*> -o <output directory> <br>

**possible “types of figure”:** <br>
csm <br>
lum-mag-veloc <br>
lum-mag-veloc-Tthresh <br>
lum-veloc-normalized <br>
lum-veloc-twostep <br>

**Comparison plots outputs:** <br>
Output files will be stored in ‘/mcmc/figures/<output directory>/<number of walkers>’. <br>
Output files will include summary figures: <br>
Comparison of fits of final models for each run-type compared, including the input observation data for that SN (all walkers, and steps after burn-in) <br>
Comparison corner plot with the posterior distributions for each run-type, overlaid on top of each-other <br>
