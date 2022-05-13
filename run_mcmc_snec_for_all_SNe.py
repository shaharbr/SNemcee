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
time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

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

def loopy_snec_mcmc(SN_name, n_steps, burn_in, n_walkers, output_dir, parameter_ranges, run_type, csm, Tthreshold, normalization, nonuniform_priors):
    # calculate range for the scaling parameter, which is derived from the uncertainty in distance
    # (in Mpc, provided in the table distances.csv')
    distances = pd.read_csv(os.path.join('SN_data', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / distance
    S_range = [1.0 - 3 * sigma_S, 1.0 + 3 * sigma_S]
    parameter_ranges['S'] = S_range
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
        flat_sampler_no_burnin = sampler_chain_flat[n_walkers * burn_in:, :]
        final_step = 2 * n_steps - 1

    mcmc_snec.chain_plots(sampler_chain, parameter_ranges, res_dir, burn_in)
    mcmc_snec.corner_plot(flat_sampler_no_burnin, parameter_ranges, res_dir)
    mcmc_snec.plot_fit_with_data(sampler_chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, SN_name, 0, Tthreshold, normalization)
    mcmc_snec.plot_fit_with_data(sampler_chain, SN_data_all, parameter_ranges,
                                       run_type, res_dir, SN_name, final_step, Tthreshold, normalization)
    np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), sampler_chain_flat, delimiter=",")

    print(sampler.chain.shape)

'''
Script can be run by calling run_mcmc_snec_for_all_SNe.py with user-provided
arguments: SN name, number of steps, number of steps and name for output directory) 
'''

def main(argv):
    SN_name = ''
    fitting_type = 'lum'
    n_steps = 500
    burn_in = 300
    n_walkers = 30
    output_dir = 'output_'+time_now
    csm = 'with'
    normalization = False
    LumThreshold = False
    nonuniform_priors = None
    arg_help = '{0} -S <SN name> -f <type of figure> -s <number of steps> -b <burn-in steps> -w <number of walkers> -o <output directory> -c <csm> -n <normalization> -Lt <luminosity threshold> -p <nonuniform priors>\n'\
               '\nargs:\n' \
               '-S SN name [required. for example: SN2017eaw]\n'\
               '-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag, combined [default: lum] \n' \
               '-s number of steps = <int> [default: 500]\n'\
               '-b burn-in steps = <int> [default: 300]\n'\
               '-w number of walkers = <int> [default: 30]\n' \
               '-o output directory name = <str> [default: output_<current time>]\n' \
               '-c csm = with, without, twostep, twostep-carryover [default: with] \n' \
               '-n normalization = True, False [default: False] \n' \
               '-Lt luminosity_threshold = True, False [default: False] \n'\
               '-p nonuniform_priors = None, <dictionary> [default: None] \n' \
               ''.format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hS:t:s:b:w:o:c:n:Lt:p", ["help", "SN=", 'fitting_type='
                                                                     "steps=", "burn_in=", "walkers=", "output_dir=",
                                                                     "csm=", "normalization=", "luminosity threshold=",
                                                                     "nonuniform_priors="])
    except:
        print(arg_help)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-S", "--SN_name"):
            SN_name = arg
        elif opt in ("-t", "--fitting_type"):
            fitting_type = arg
        elif opt in ("-s", "--steps"):
            n_steps = int(arg)
        elif opt in ("-b", "--burn_in"):
            burn_in = int(arg)
        elif opt in ("-w", "--walkers"):
            n_walkers = int(arg)
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-c", "--csm"):
            csm = arg
        elif opt in ("-n", "--normalization"):
            normalization = arg
        elif opt in ("-Lt", "--luminosity_threshold"):
            LumThreshold = arg
        elif opt in ("-p", "--nonuniform_priors"):
            nonuniform_priors = arg

    print('SN name:', SN_name)
    print('fitting type:', fitting_type)
    print('steps:', n_steps)
    print('walkers:', n_walkers)
    print('burn in:', burn_in)
    print('output directory:', output_dir)
    print('csm:', csm)
    print('normalization:', normalization)
    print('luminosity threshold:', LumThreshold)
    print('nonuniform priors:', nonuniform_priors)

    loopy_snec_mcmc(SN_name, n_steps, burn_in, n_walkers, output_dir,
                    parameter_ranges, fitting_type, csm,
                    {'lum': LumThreshold, 'mag': True, 'veloc': True},
                    normalization,
                    nonuniform_priors)



if __name__ == "__main__":
    main(sys.argv)

