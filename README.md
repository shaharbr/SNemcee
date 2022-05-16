# SNemcee

**Dependencies:** <br>
Python libraries: <br>
Corner <br>
Emcee <br>
Pandas <br>
Numpy <br>
Datetime <br>
scipy <br>


**Input supernova data:** <br>
Supernova data files (observations) should be stored in the directory ‘SN_data’, with the following names: <br>
<SN_name>_veloc <br>
<SN_name>_mag <br>
<SN_name>_lum <br>
(See example files for expected format) <br>

In addition, there should be a table in ‘SN_data’ that contains the distance (in Mpc) <br>


**SNEC model grid:** <br>
Should be downloaded from:
https://drive.google.com/file/d/1-7ki0ADkq2_Y8ffyzWMT42_VuvXQK69p/view?usp=sharing

And then added to the main directory and extracted using
tar -xf SNEC_models.tar.gz
The result should be the directory SNemcee/SNEC_models which contains 2,469 sub-directories with names like 'M9.0_Ni0.001_E0.1_Mix2.0_R1000_K10'

**Running SNemcee:** <br>
To run the SNemcee program on your SN of choice, use command: <br>

python3 run_mcmc_snec_for_all_SNe.py -S <SN name> -t <fitting_type> -s <number of steps> -b <burn-in steps> -w <number of walkers> -o <output directory> -c <csm> -n <normalization> -l <luminosity threshold> -p <nonuniform priors> <br>

-S SN name [required. for example: SN2017eaw]<br>
-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag, combined [default: lum]<br>
-s number of steps = <int> [default: 500]<br>
-b burn-in steps = <int> [default: 300]<br>
-w number of walkers = <int> [default: 30]<br>
-o <output directory> = <string> [default: output_<current time>]<br>
-c csm = with, without, twostep, twostep-carryover [default: with]<br>
-n normalization = True, False [default: False]<br>
-l luminosity_threshold = True, False [default: False]<br>
-p nonuniform_priors = None, <dictionary> [default: None]<br>


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

python3 plot_snec_fits_wrapper.py -S <SN name> -s <number of steps> -f <figure_type> -o <output directory> -t <fitting type> -c <csm> -n <normalization> -l <luminosity threshold> <br>


-S SN name [required. for example: SN2017eaw]<br>
-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag, combined [default: lum]<br>
-s number of steps = <int> [default: 500]<br>
-f figure type = single fit-type plots: lum/veloc/mag [or any combination of those, separated by -] OR corner OR comparison plots: csm_comparison/lum-mag-veloc_comparison/lum-mag-veloc-Tthresh_comparison/lum-veloc-normalized/lum-veloc-twostep [required]<br>
-o <output directory> = <string> [default: output_<current time>]<br>
-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag, combined [default: False] <br>
-c csm = with, without, twostep, twostep-carryover [default: with]<br>
-n normalization = True, False [default: False]<br>
-l luminosity_threshold = True, False [default: False]<br>


**Comparison plots outputs:** <br>
Output files will be stored in ‘/mcmc/figures/<output directory>/<number of walkers>’. <br>
Output files will include summary figures: <br>
Comparison of fits of final models for each run-type compared, including the input observation data for that SN (all walkers, and steps after burn-in) <br>
Comparison corner plot with the posterior distributions for each run-type, overlaid on top of each-other <br>
    