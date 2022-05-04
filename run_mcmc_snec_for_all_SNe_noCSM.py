import numpy as np
import pandas as pd
from pathlib import Path
import os
import mcmc_snec_noCSM as mcmc_snec
import datetime

date = '2022_02_23'
n_steps = 300


def loopy_snec_mcmc(SN_name, parameter_ranges, run_type):
    distances = pd.read_csv(os.path.join('results', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / distance
    S_range = [1.0 - sigma_S, 1.0 + sigma_S]
    parameter_ranges['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
    time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

    n_walkers = 30
    n_params = 6
    burn_in = int(n_steps * 0.4)
    res_dir = os.path.join('mcmc_results', date, str(n_steps) + 'step',
                           run_type + '_noCSM_' + SN_name)

    Path(res_dir).mkdir(parents=True, exist_ok=True)
    mcmc_snec.write_params_file(parameter_ranges, S_range, SN_name, n_walkers, n_steps, n_params, burn_in, time_now,
                                res_dir)
    SN_data_all = {}
    if 'lum' in run_type:
        # import SN bolometric lum data
        data_lum = mcmc_snec.import_lum(SN_name)
        SN_data_all['lum'] = data_lum
    if 'veloc' in run_type:
        # import SN photospheric velocities data
        data_veloc = mcmc_snec.import_veloc(SN_name)
        SN_data_all['veloc'] = data_veloc
    if 'mag' in run_type:
        # import SN mag data
        data_mag = mcmc_snec.import_mag(SN_name)
        SN_data_all['mag'] = data_mag

    '''
    # running code #
    '''

    sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                         run_type)

    flat_sampler = sampler.get_chain(flat=True)
    np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

    flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
    np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")

    mcmc_snec.chain_plots(sampler.chain, parameter_ranges, res_dir, burn_in)
    mcmc_snec.corner_plot(sampler.get_chain(flat=True)[n_walkers * burn_in:, :], parameter_ranges, res_dir)

    mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, n_walkers, SN_name, n_steps - 1)
    mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, n_walkers, SN_name, 0)

    print(sampler.chain.shape)


Mzams_range = [9.0, 10.0, 11.0, 13.0, 15.0, 17.0]
Ni_range = [0.001, 0.02, 0.07, 0.12]
E_final_range = [0.1, 0.3, 0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
T_range = [-5, 2]
parameter_ranges = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                    'Mix': Mix_range, 'S': [], 'T': T_range}

SN_list_lum = ['SN2004a', 'SN2004et_d5.9', 'SN2012aw', 'SN2012ec', 'SN2017eaw', 'SN2018aoq', 'SN2005cs', 'SN2008bk']
SN_list_veloc = ['SN2012aw', 'SN2017eaw', 'SN2005cs', 'SN2008bk']  # veloc
SN_list_mag = ['SN2017eaw', 'SN2004et_d5.9']  # mag
SN_list = ['SN2020bij']

# for run_type in ['lum-veloc']:
#     for SN_name in ['SN2017eaw', 'SN2005cs']:
#         loopy_snec_mcmc(SN_name, parameter_ranges, run_type)
# for run_type in ['veloc', 'mag', 'lum']:
#   for SN_name in ['SN2017eaw']:
#      loopy_snec_mcmc(SN_name, parameter_ranges, run_type)

for run_type in ['lum-veloc-normalized']:
    for SN_name in SN_list:
        loopy_snec_mcmc(SN_name, parameter_ranges, run_type)
for run_type in ['mag-veloc-normalized']:
    for SN_name in SN_list:
        loopy_snec_mcmc(SN_name, parameter_ranges, run_type)