import numpy as np
import pandas as pd
from pathlib import Path
import os
import mcmc_snec
import mcmc_snec_noCSM
import datetime

date = '2022_02_23'
n_steps = 50

def loopy_snec_mcmc(SN_name, parameter_ranges, run_type, priors_carryover=None):
    distances = pd.read_csv(os.path.join('SN_data', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / distance
    S_range = [1.0 - sigma_S, 1.0 + sigma_S]
    parameter_ranges['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
    parameter_ranges_no_CSM = parameter_ranges.copy()
    parameter_ranges_no_CSM.pop('R')
    parameter_ranges_no_CSM.pop('K')
    time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

    n_walkers = 16
    n_params = 8
    burn_in = int(n_steps * 0.4)

    res_dir = os.path.join('mcmc_results', date, str(n_steps) + 'step',
                           run_type + '_twostep_priors' +str(priors_carryover) + SN_name)

    Path(res_dir).mkdir(parents=True, exist_ok=True)
    mcmc_snec.write_params_file(parameter_ranges, S_range, SN_name, n_walkers, n_steps, n_params, burn_in, time_now,
                                res_dir)
    SN_data_all = {}
    if 'lum' in run_type:
        # import SN bolometric lum SNEC_models
        data_lum = mcmc_snec.import_lum(SN_name)
        SN_data_all['lum'] = data_lum
    if 'veloc' in run_type:
        # import SN photospheric velocities SNEC_models
        data_veloc = mcmc_snec.import_veloc(SN_name)
        SN_data_all['veloc'] = data_veloc
    if 'mag' in run_type:
        # import SN mag SNEC_models
        data_mag = mcmc_snec.import_mag(SN_name)
        SN_data_all['mag'] = data_mag

    '''
    # running code #
    '''

    """
    first step: fit to Mzams, Ni, E, Mix, S, T, without CSM. With gaussian prior for S, uniform for all others.
    """
    sampler_first_step = mcmc_snec_noCSM.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges_no_CSM,
                                                          run_type)

    # adding zero CSM columns to the sampler chain of the first step
    sampler_first_step_chain_flat = sampler_first_step.get_chain(flat=True)
    sampler_first_step_chain_flat_noburnin = sampler_first_step.get_chain(discard=burn_in, flat=True)
    sampler_first_step_chain = sampler_first_step.chain
    sampler_first_step_chain_flat_CSM = np.zeros((n_steps * n_walkers, 8))
    sampler_first_step_chain_flat_CSM_noburnin = np.zeros(((n_steps - burn_in) * n_walkers, 8))
    sampler_first_step_chain_CSM = np.zeros((n_walkers, n_steps, 8))
    for i in range(6):
        if i < 3:
            sampler_first_step_chain_flat_CSM[:, i] = sampler_first_step_chain_flat[:, i]
            sampler_first_step_chain_flat_CSM_noburnin[:, i] = sampler_first_step_chain_flat_noburnin[:, i]
            sampler_first_step_chain_CSM[:, :, i] = sampler_first_step_chain[:, :, i]
        else:
            sampler_first_step_chain_flat_CSM[:, i + 2] = sampler_first_step_chain_flat[:, i]
            sampler_first_step_chain_flat_CSM_noburnin[:, i + 2] = sampler_first_step_chain_flat_noburnin[:, i]
            sampler_first_step_chain_CSM[:, :, i + 2] = sampler_first_step_chain[:, :, i]

    """
    second step: fit to K and R, with priors according to SN_data from step one.
    """

    flat_sampler_no_burnin = pd.DataFrame(sampler_first_step_chain_flat_CSM[n_walkers * burn_in:, :])
    flat_sampler_no_burnin.columns = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']

    if priors_carryover is not None:
        nonuniform_priors = {}
        for param in ['Mzams', 'Ni', 'E', 'Mix', 'S', 'T']:
            nonuniform_priors[param] = {
                'polynomial': mcmc_snec.polyfit_to_distribution(flat_sampler_no_burnin[param], res_dir)}
    else:
        nonuniform_priors = None

    last_walkers_step_one = sampler_first_step_chain_flat_CSM[-n_walkers:, :]
    # change R and K to initial guesses of uniform dist within range
    last_walkers_step_one[:, 3] = np.random.rand(n_walkers) * \
                                  (parameter_ranges['R'][-1] - parameter_ranges['R'][0]) + parameter_ranges['R'][0]
    last_walkers_step_one[:, 4] = np.random.rand(n_walkers) * \
                                  (parameter_ranges['K'][-1] - parameter_ranges['K'][0]) + parameter_ranges['K'][0]

    sampler_second_step = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                                     run_type,
                                                     nonuniform_priors, init_guesses=last_walkers_step_one)

    sampler_chain = np.concatenate((sampler_first_step_chain_CSM, sampler_second_step.chain), axis=1)

    sampler_chain_flat = np.concatenate((
        sampler_first_step_chain_flat_CSM,
        sampler_second_step.get_chain(flat=True)), axis=0)

    sampler_chain_flat_noburnin = np.concatenate((
        sampler_first_step_chain_flat_CSM_noburnin,
        sampler_second_step.get_chain(flat=True)), axis=0)

    np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), sampler_chain_flat, delimiter=",")
    np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), sampler_chain_flat_noburnin, delimiter=",")

    mcmc_snec.chain_plots(sampler_chain, parameter_ranges, res_dir, burn_in, first_stage_steps=n_steps)
    mcmc_snec.corner_plot(sampler_chain_flat[n_walkers * (n_steps + burn_in):, :], parameter_ranges, res_dir)
    mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, n_walkers, SN_name, n_steps * 2 - 1)
    print(sampler_chain.shape)


Mzams_range = [9.0, 10.0, 11.0, 13.0, 15.0, 17.0]
Ni_range = [0.001, 0.02, 0.07, 0.12]
E_final_range = [0.1, 0.3, 0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
R_range = [0, 500, 1000]
K_range = [0, 10, 30, 60]
T_range = [-5, 2]
parameter_ranges = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                    'R': R_range, 'K': K_range, 'Mix': Mix_range, 'S': [], 'T': T_range}

SN_list_lum = ['SN2004a', 'SN2004et_d5.9', 'SN2012aw', 'SN2012ec', 'SN2017eaw', 'SN2018aoq', 'SN2005cs', 'SN2008bk']
SN_list_veloc = ['SN2012aw', 'SN2017eaw','SN2005cs', 'SN2008bk']  # veloc
SN_list_mag = ['SN2017eaw', 'SN2004et_d5.9']  # mag
SN_list = ['SN2017eaw']

for run_type in ['veloc', 'lum-veloc', 'lum-veloc-normalized']:
    for SN_name in SN_list:
        for priors_carryover in [None, True]:
            loopy_snec_mcmc(SN_name, parameter_ranges, run_type, priors_carryover)
for run_type in ['mag-veloc', 'mag-veloc-normalized']:
    for SN_name in SN_list:
        for priors_carryover in [None, True]:
            loopy_snec_mcmc(SN_name, parameter_ranges, run_type, priors_carryover)
