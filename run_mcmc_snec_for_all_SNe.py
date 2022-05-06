'''

a
This script runs the mcmc_snec code based on the user-provided arguments:
Supernova (SN) name
number  of steps
number of walkers
output directory (usually would be the date)
by default, number of burn-in steps is 40% of total step number
'''

import numpy as np
import pandas as pd
import os
import mcmc_snec
import datetime
import sys, getopt


'''
Lists specifying the range of parameters to allow the mcmc to explore,
based on the extent of the SNEC grid available
'''
Mzams_range = [9.0, 10.0, 11.0, 13.0, 15.0, 17.0]  # Progenitor ZAMS mass (solar mass)
Ni_range = [0.001, 0.02, 0.07, 0.12]  # Ni mass (solar mass)
E_final_range = [0.1, 0.3, 0.5, 0.9, 1.3, 1.7]  # Explosion energy (foe=10^51 erg)
Mix_range = [2.0, 8.0]  # Ni mixing (solar mass)
R_range = [0, 500, 1000]
K_range = [0, 10, 30, 60]
# R_range = [0,0]
# K_range = [0,0]
T_range = [-5, 2]  # time shift from assumed day of explosion (days)
parameter_ranges = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                    'R': R_range, 'K': K_range, 'Mix': Mix_range, 'S': [], 'T': T_range}


'''
The function loopy_snec_mcmc imports the SNEC_models for the SN, and then runs
mcmc_snec.py with the user-provided arguments: SN, step number, walker number 
as well as the run-type arguments determined by the main funciton:
*run_type: 'lum', 'veloc', 'mag', 'lum-veloc', 'mag-veloc', 'lum-mag', 'combined'
[which observations are used to calculate the fitting score in the likelihood function]
*csm: with, without, twostep, twostep_carryover
*Tthreshold: dictionary {'lum': True/False, 'mag': True/False, 'veloc': True/False}
*normalization: True or False
*nonuniform_priors: dictionary with keys being parameter names, and their values 
being more dictionaries with gaussian or polynomial distribution parameters. 
For example:
{'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}
{'E': {'polynomial': np.poly1d(np.polyfit()) output}

The function outputs the MCMC SN_data, figures and running parameters used to
the output directory. Output files include:
1) Corner plot of posterior distributions for each parameter (all walkers, and steps after burn-in)
2) Fits of final models to the input observation SNEC_models for that SN (all walkers, and steps after burn-in)
3) Chain plots for each parameter
4) Flat_sampler.csv: MCMC sampler output table - a 2D table, where columns are the parameters, and rows are walkers and steps (by the hierarchy of: steps, and within them walkers)
5) final_results.csv: Summary statistics for the final posterior distributions for each parameter: mean, 16% percentile and 85% percentile.
6) run_parameters.csv: Record of the input arguments used in that run



'''

def loopy_snec_mcmc(SN_name, n_steps, n_walkers, output_dir, parameter_ranges, run_type, csm, Tthreshold, normalization, nonuniform_priors):
    # default burn_in steps are 40% of total numebr of steps
    burn_in = int(n_steps * 0.4)
    # calculate range for the scaling parameter, which is derived from the uncertainty in distance
    # (in Mpc, provided in the table distances.csv')
    distances = pd.read_csv(os.path.join('SN_data', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / distance
    S_range = [1.0 - sigma_S, 1.0 + sigma_S]
    parameter_ranges['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
    time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = SN_name + '_' +  run_type + '_csm-' + str(csm) + '_' + 'normalized' + str(normalization)
    # TODO think about how user should provide the nonuniform_priors
    if nonuniform_priors is not None:
        nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

    SN_data_all = {}
    if 'lum' in run_type or 'combined' in run_type:
        # import SN bolometric lum SNEC_models
        data_lum = mcmc_snec.import_lum(SN_name)
        SN_data_all['lum'] = data_lum
        filename += '_' + 'TreshLum' + str(Tthreshold['lum'])
    if 'mag' in run_type or 'combined' in run_type:
        # import SN mag SNEC_models
        data_mag = mcmc_snec.import_mag(SN_name)
        SN_data_all['mag'] = data_mag
        filename += '_' + 'TreshMag' + str(Tthreshold['mag'])
    if 'veloc' in run_type or 'combined' in run_type:
        # import SN photospheric velocities SNEC_models
        data_veloc = mcmc_snec.import_veloc(SN_name)
        SN_data_all['veloc'] = data_veloc
        filename += '_' + 'TreshVeloc' + str(Tthreshold['veloc'])

    res_dir = os.path.join('mcmc_results', output_dir, str(n_steps) + 'step', filename)
    if not os.path.exists(res_dir):
        if not os.path.exists(os.path.join('mcmc_results', output_dir, str(n_steps) + 'step')):
            if not os.path.exists(os.path.join('mcmc_results', output_dir)):
                os.mkdir(os.path.join('mcmc_results', output_dir))
            os.mkdir(os.path.join('mcmc_results', output_dir, str(n_steps) + 'step'))
        os.mkdir(res_dir)
    mcmc_snec.write_params_file(parameter_ranges, S_range, SN_name, n_walkers, n_steps,
                                csm, Tthreshold, normalization, burn_in, time_now, res_dir)
    '''
    # running code #
    '''
    final_step = n_steps - 1
    if csm == 'with':
        sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                             run_type, True, Tthreshold, normalization, nonuniform_priors)
        sampler_chain = sampler.chain
    elif csm == 'without' or csm == 'twostep' or csm == 'twostep-carryover':
        sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                             run_type, False, Tthreshold, normalization, nonuniform_priors)
        sampler_chain = sampler.chain
        sampler_chain = np.insert(sampler_chain, 3, np.zeros((2, 1, 1)), axis=2)

    # make flat (2D) chain matrix
    s = list(sampler_chain.shape[1:])
    s[0] = np.prod(sampler_chain.shape[:2])
    sampler_chain_flat = sampler_chain.reshape(s)
    # make flat chain without the burn-in steps
    flat_sampler_no_burnin = sampler_chain_flat[n_walkers * burn_in:, :]

    if csm == 'twostep' or csm == 'twostep-carryover':
        # make a flat chain matrix as dataframe with column headers
        flat_sampler_no_burnin_df = pd.DataFrame(flat_sampler_no_burnin)
        flat_sampler_no_burnin_df.columns = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']

        # get the positions of the walkers in the last step of the no-CSM stage
        last_walkers_noCSM = flat_sampler_no_burnin[-n_walkers:, :]
        # change the R and K of that last step to initial guesses of uniform distribution within their range
        last_walkers_noCSM[:, 3] = np.random.rand(n_walkers) * \
                                      (parameter_ranges['R'][-1] - parameter_ranges['R'][0]) + \
                                      parameter_ranges['R'][0]
        last_walkers_noCSM[:, 4] = np.random.rand(n_walkers) * \
                                      (parameter_ranges['K'][-1] - parameter_ranges['K'][0]) + \
                                      parameter_ranges['K'][0]
        # if we're doing carryover of the non-CSM results as priors to the second stage, extract their ditributions
        # in the last step and set as priors
        if csm == 'twostep':
            nonuniform_priors = None
        elif csm == 'twostep-carryover':
            nonuniform_priors = {}
            for param in ['Mzams', 'Ni', 'E', 'Mix', 'S', 'T']:
                nonuniform_priors[param] = {
                    'polynomial': mcmc_snec.polyfit_to_distribution(flat_sampler_no_burnin_df[param], res_dir)}
        # run the mcmc again, with walkers starting from their positions at the end of the first stage, and with or
        # without carryover as priors
        sampler_second_step = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                                         run_type, True, Tthreshold, normalization, nonuniform_priors,
                                                         init_guesses=last_walkers_noCSM)
        # concatenate the second stage sampler to the first one
        sampler_chain = np.concatenate((sampler_chain, sampler_second_step.chain), axis=1)
        sampler_chain_flat = np.concatenate((sampler_chain_flat, sampler_second_step.get_chain(flat=True)), axis=0)
        final_step = 2 * n_steps - 1

    # TODO change name of function - its for any plots (not only lum)
    mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, SN_name, 0, Tthreshold, normalization)
    mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, SN_name, final_step, Tthreshold, normalization)
    np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), sampler_chain_flat, delimiter=",")

    print(sampler.chain.shape)

'''
Script can be run by calling run_mcmc_snec_for_all_SNe.py with user-provided
arguments: SN name, number of steps, number of steps and name for output directory) 
'''

def main(argv):
    SN_name = ""
    n_steps = 0
    n_walkers = 0
    output_dir = ""
    arg_help = "{0} -S <SN name> -s <number of steps> -w <number of walkers> -o <output directory>".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hS:s:w:o:", ["help", "SN=",
                                                         "steps=", "walkers=", "output_dir="])
    except:
        print(arg_help)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-S", "--SN_name"):
            SN_name = arg
        elif opt in ("-s", "--steps"):
            n_steps = int(arg)
        elif opt in ("-w", "--walkers"):
            n_walkers = int(arg)
        elif opt in ("-o", "--output_dir"):
            output_dir = arg

    print('SN name:', SN_name)
    print('steps:', n_steps)
    print('walkers:', n_walkers)
    print('output directory:', output_dir)
    # if os.path.exists(os.path.join('SN_data', SN_name + '_veloc')):
    #     for normalization in [True, False]:
    #         for LumT in [False, True]:
    #             loopy_snec_mcmc(SN_name, n_steps, n_walkers, output_dir,
    #                             parameter_ranges, 'lum-veloc',
    #                             csm='with',
    #                             Tthreshold={'lum': LumT, 'mag': True, 'veloc': True},
    #                             normalization=normalization,
    #                             nonuniform_priors=None)
    #     if os.path.exists(os.path.join('SN_data', SN_name + '_mag')):
    #         for normalization in [True, False]:
    #             for LumT in [False, True]:
    #                 loopy_snec_mcmc(SN_name, n_steps, n_walkers, output_dir,
    #                                 parameter_ranges, 'mag-veloc',
    #                                 csm='with',
    #                                 Tthreshold={'lum': LumT, 'mag': True, 'veloc': True},
    #                                 normalization=normalization,
    #                                 nonuniform_priors=None)
    for csm in ['twostep', 'twostep-carryover', 'with', 'without']:
        loopy_snec_mcmc(SN_name, n_steps, n_walkers, output_dir,
                        parameter_ranges, 'lum',
                        csm,
                        Tthreshold={'lum': False, 'mag': True, 'veloc': True},
                        normalization=False,
                        nonuniform_priors=None)


if __name__ == "__main__":
    main(sys.argv)


# user input parameters:
# output_dir = '2022_03_23'
# SN_list = ['SN2020bij']
# n_steps = 500
# n_walkers = 50


# SN_list_lum = ['SN2004a', 'SN2004et_d5.9', 'SN2012aw', 'SN2012ec', 'SN2017eaw','SN2018aoq', 'SN2005cs', 'SN2008bk']
# SN_list_veloc = ['SN2012aw', 'SN2004et_d5.9', 'SN2017eaw','SN2005cs', 'SN2008bk']  # veloc
# SN_list_mag = ['SN2017eaw', 'SN2004et_d5.9']  # mag


# Tthreshold = {'lum': False, 'veloc': True, 'mag': True}
# normalization = False
# nonuniform_priors = None
# csm = True


# for csm in [True, False]:
#     loopy_snec_mcmc(SN_name, parameter_ranges, 'lum',
#                     csm,
#                     Tthreshold={'lum': False, 'mag': True, 'veloc': True},
#                     normalization=False,
#                     nonuniform_priors=None)
# for run_type in ['lum-veloc', 'mag-veloc', 'combined']:
#    for normalization in [True, False]:
#        for LumT in [False, True]:
#            loopy_snec_mcmc(SN_name, parameter_ranges, run_type,
#                            csm=True,
#                            Tthreshold={'lum': LumT, 'mag': True, 'veloc': True},
#                            normalization=normalization,
#                            nonuniform_priors=None)


#for run_type in ['veloc', 'lum-veloc', 'lum-veloc-normalized']:
 #   for SN_name in SN_list_veloc:
  #      loopy_snec_mcmc(SN_name, parameter_ranges, run_type)
#for run_type in ['mag-veloc', 'mag-veloc-normalized']:
 #   for SN_name in SN_list_mag:
  #      loopy_snec_mcmc(SN_name, parameter_ranges, run_type)
