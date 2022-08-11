import numpy as np
import emcee
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import snec_model_interpolator as interp
import pandas as pd
import os
import csv
from scipy.stats import chi2
import copy
from numpy import trapz
import models as mod
import datetime


'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'

colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}

# multiplicative factor for SNEC's predicted photospheric velocities to Fe II velocities,
# found to be equal about 1.4, by fitting to all the SNe in our sample
T_thresh = 10 ** 3.75
models = {}


def polyfit_to_distribution(array_walker_results, res_dir):
    name = array_walker_results.name
    counts, bin_edges = np.histogram(array_walker_results, bins=20)
    fig, ax = plt.subplots()
    ax.hist(array_walker_results, density=True, bins=20)
    bin_widths = np.diff(bin_edges)
    x = bin_edges[:-1] + (bin_widths / 2)
    y = counts
    area = trapz(y, dx=bin_widths[0])
    y = y / area
    polymod = np.poly1d(np.polyfit(x, y, deg=12))
    ax.plot(x, y, color='black')
    dense_x = np.arange(np.min(x), np.max(x), (np.max(x)-np.min(x))/50)
    ax.plot(dense_x, [polymod(i) for i in dense_x], color='orange')
    fig.savefig(os.path.join(res_dir, str(name)+'_first_step_dist.png'))
    return polymod


def load_surrounding_models(requested, ranges_dict, fitting_type, LumTthreshold=False):
    mod.load_surrounding_models_local(models, requested, ranges_dict, fitting_type, LumTthreshold)


def theta_in_range(theta, ranges_dict):
    Mzams = theta[0]
    E = theta[2]
    if Mzams > 11 and E < 0.5:
        return False
    else:
        ranges_list = dict_to_list(ranges_dict)
        truth_value = True
        for i in range(len(ranges_list)):
            truth_value = truth_value and\
                          (ranges_list[i][0] <= theta[i] <= ranges_list[i][-1])
        return truth_value


def log_prior(theta, ranges_dict, nonuniform_priors=None):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form of (example):
    ranges_dict =  {Mzams : [13.0, 14.0, 15.0, 16.0],
                    Ni : [0.02, 0.07, 0.12, 0.17],
                    E : [0.7, 0.9, 1.1, 1.3, 1.5],
                    Mix : [5.0, 8.0],
                    R : [0.1, 0.1],
                    K : [0.1, 0.1],
                    S : [0.8, 1.2],
                    T : [-15, +15]}
    """
    if not theta_in_range(theta, ranges_dict):
        return -np.inf
    else:
        if nonuniform_priors is None:
            # flat probability for all means p=1 for all, and log_p=0, so sum is also 0.
            return 0.0
        else:
            param_dict = {'Mzams': theta[0], 'Ni': theta[1], 'E': theta[2], 'R': theta[3],
                          'K': theta[4], 'Mix': theta[5], 'S': theta[6], 'T': theta[7]}
            log_prior = 0.0
            for param in list(nonuniform_priors.keys()):
                if nonuniform_priors[param] == 'gaussian':
                    mu = nonuniform_priors[param]['gaussian']['mu']
                    sigma = nonuniform_priors[param]['gaussian']['sigma']
                    log_prior += np.log(1.0 / (np.sqrt(2 * np.pi) * sigma)) - 0.5 * (param_dict[param] - mu) ** 2 / sigma ** 2
                    # print('logprior', log_prior)
                if nonuniform_priors[param] == 'polynomial':
                    log_prior += np.log(nonuniform_priors[param]['polynomial'](param_dict[param]))
            return log_prior


def dict_to_list(dict):
    l = []
    for key in dict.keys():
        l.append(dict[key])
    return l


def chi_square_norm(x, y, dy, x_fit, y_fit):
    y_fit_on_data_x = np.interp(x, x_fit, y_fit)
    return np.sum(((y - y_fit_on_data_x) / dy) ** 2)


def interp_yfit(theta, ranges_dict, fitting_type, data_x, filter=False):
    surrounding_values = mod.get_surrouding_values(theta[0:6], ranges_dict)
    load_surrounding_models(theta[0:6], ranges_dict, fitting_type)
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models[fitting_type], data_x, filter)
    return y_fit


# TODO this step can introduce some bias - SNe that don't have early velocities
#  will select against models that cool fast
def temp_thresh_cutoff(requested_theta, ranges_dict, data_x):
    temp_fit = interp_yfit(requested_theta, ranges_dict, 'temp', data_x)
    max_x = 0
    if not isinstance(temp_fit, str):
        if len(temp_fit[temp_fit >= T_thresh]) > 0:
            min_temp_above_Tthresh = np.min(temp_fit[temp_fit >= T_thresh])
            x_out = data_x[temp_fit > min_temp_above_Tthresh]
            max_x = np.max(x_out)
    return max_x, temp_fit


def calc_likelihood(data_x, data_y, data_dy, y_fit, normalization):
    # calculate the log likelihood
    df = len(data_y) - 1
    if normalization:
        norm = len(data_y)
    else:
        norm = 1
    return chi2.logpdf(chi_square_norm(data_x, data_y, data_dy, data_x, y_fit), df) / norm


def calc_lum_likelihood(theta, data_dict, ranges_dict, normalization, LumTthreshold=False):
    data = data_dict['lum']
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    if LumTthreshold:
        max_x, temp_fit = temp_thresh_cutoff(theta[0:6], ranges_dict, data_x_moved)
        data = data.loc[data_x_moved <= max_x]
        data_x_moved = data_x_moved[data_x_moved <= max_x]
    data_y = data['Lum']
    data_dy = data['dLum0']
    y_fit = interp_yfit(theta, ranges_dict, 'lum', data_x_moved)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        likeli = calc_likelihood(data_x_moved, data_y, data_dy, y_fit, normalization)
        return likeli
    else:
        # print('impossible SN')
        return - np.inf


def calc_veloc_likelihood(theta, data_dict, ranges_dict, normalization):
    data = data_dict['veloc']
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    # apply temperature threshold on time
    max_x, temp_fit = temp_thresh_cutoff(theta[0:6], ranges_dict, data_x_moved)
    data = data.loc[data_x_moved <= max_x]
    data_x_moved = data_x_moved[data_x_moved <= max_x]
    # veloc data
    data_y = data['veloc']
    data_dy = data['dveloc']
    y_fit = interp_yfit(theta, ranges_dict, 'veloc', data_x_moved)
    if not isinstance(y_fit, str):
        # calculate the log likelihood
        likeli = calc_likelihood(data_x_moved, data_y, data_dy, y_fit, normalization)
        return likeli
    else:
        # print('impossible SN')
        return - np.inf


def calc_mag_likelihood(theta, data_dict, ranges_dict, normalization):
    data = data_dict['mag']
    likeli = 0
    any_filter_data = False
    y_fit = {}
    filters = list(data['filter'].unique())
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    # apply temperature threshold on time
    max_x, temp_fit = temp_thresh_cutoff(theta[0:6], ranges_dict, data_x_moved)
    data = data.loc[data_x_moved <= max_x]
    # mag data
    for filt in filters:
        data_filt = data.loc[data['filter'] == filt]
        data_x_filt_moved = data_filt['t_from_discovery'] - theta[7]
        data_y_filt = data_filt['abs_mag']
        data_dy_filt = data_filt['dmag']
        y_fit[filt] = interp_yfit(theta, ranges_dict, 'mag', data_x_filt_moved, filter=filt)
        if not isinstance(y_fit[filt], str):
            # multiply whole graph by scaling factor
            y_fit[filt] = y_fit[filt] -2.5*np.log10(theta[6])
            # calculate the log likelihood
            likeli += calc_likelihood(data_x_filt_moved, data_y_filt, data_dy_filt, y_fit[filt],
                                      normalization)
            any_filter_data = True
    if any_filter_data:
        return likeli
    else:
        return - np.inf


def log_likelihood(theta, data, ranges_dict, fitting_type, LumTthreshold, normalization):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form described above

    fitting_type : lum, mag, veloc
    """
    # print(theta)
    # ranges_list = dict_to_list(ranges_dict)
    # print(ranges_list)
    if theta_in_range(theta, ranges_dict):
        # print('ok SN')
        load_surrounding_models(theta[0:6], ranges_dict, fitting_type, LumTthreshold)
        args = [theta, data, ranges_dict, normalization]
        # TODO make this run any combination of fitting type based on "if in", for, add
        if fitting_type == 'lum':
            args.append(LumTthreshold)
            log_likeli = calc_lum_likelihood(*args)
        elif fitting_type == 'mag':
            log_likeli = calc_mag_likelihood(*args)
        elif fitting_type == 'veloc':
            log_likeli = calc_veloc_likelihood(*args)
        elif fitting_type == 'lum-veloc':
            log_likeli = calc_lum_likelihood(*args) + \
                         calc_veloc_likelihood(*args)
        elif fitting_type == 'mag-veloc':
            log_likeli = calc_mag_likelihood(*args) + \
                         calc_veloc_likelihood(*args)
        elif fitting_type == 'combined':
            log_likeli = calc_lum_likelihood(*args) + \
                         calc_mag_likelihood(*args) + \
                         calc_veloc_likelihood(*args)
        else:
            print('fitting_type should be: lum, mag, veloc, lum-veloc, lum-veloc_normalized, mag-veloc, mag-veloc-normalized, combined or combined-normalized')
    else:
        # print('log lik')
        # print(theta)
        # print(ranges_dict)
        # print('out of range')
        return - np.inf  # just a very big number so it won't go past the edge values
    # print('loglik', log_likeli)
    return log_likeli


def log_posterior(theta, data, ranges_dict, fitting_type, csm, LumTthreshold, normalization, nonuniform_priors):
    if not csm:
        theta =list(np.insert(theta, 3, np.zeros((2,)), axis=0))
    if theta_in_range(theta, ranges_dict):
        lp = log_prior(theta, ranges_dict, nonuniform_priors)
        ll = log_likelihood(theta, data, ranges_dict, fitting_type, LumTthreshold, normalization)
        log_post = lp + ll
        # print('logpost', log_post)
        return log_post
    else:
        # print('log pos')
        # print(theta)
        # print(ranges_dict)
        # print('out of range')
        return - np.inf  # just a very big number so it won't go past the edge values


def initial_guesses(ranges_dict, n_walkers, csm):
    guesses = []
    if csm:
        for param in list(ranges_dict.keys()):
            guesses.append(np.random.rand(n_walkers) *
                           (ranges_dict[param][-1] - ranges_dict[param][0]) + ranges_dict[param][0])
    else:
        for param in list(ranges_dict.keys()):
            if param != 'R' and param != 'K':
                guesses.append(np.random.rand(n_walkers) *
                               (ranges_dict[param][-1] - ranges_dict[param][0]) + ranges_dict[param][0])
    return np.array(guesses).T


def load_SN_data(run_type, SN_name):
    SN_data_all = {}
    print(run_type)
    if 'lum' in run_type:
        # import SN bolometric lum SNEC_models
        data_lum = import_lum(SN_name)
        SN_data_all['lum'] = data_lum
    if 'mag' in run_type:
        # import SN mag SNEC_models
        data_mag = import_mag(SN_name)
        SN_data_all['mag'] = data_mag
    if 'veloc' in run_type:
        # import SN photospheric velocities SNEC_models
        data_veloc = import_veloc(SN_name)
        SN_data_all['veloc'] = data_veloc
    print(SN_data_all)
    return SN_data_all


def S_prior(SN_name):
    # calculate range for the scaling parameter, which is derived from the uncertainty in distance
    # (in Mpc, provided in the table distances.csv')
    distances = pd.read_csv(os.path.join('SN_data', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / distance
    S_range = [1.0 - 3 * sigma_S, 1.0 + 3 * sigma_S]
    return S_range, sigma_S


def save_param_results(sampler_chain, ranges_dict, step, res_dir):
    params = list(ranges_dict.keys())
    dict = {}
    for i, param in enumerate(params):
        last_results = sampler_chain[:, step:, i]
        avg = np.mean(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[param] = avg
        dict[param +'_lower'] = avg - sigma_lower
        dict[param + '_upper'] = sigma_upper - avg
    with open(os.path.join(res_dir, 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)
    return dict


def write_params_file(parameter_ranges, SN_name, n_walkers, n_steps, csm, LumTthreshold, normalization, burn_in, time_now, output_dir):
    run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                                           'Mzams_range': str(parameter_ranges['Mzams']),
                                           'Ni_range': str(parameter_ranges['Ni']),
                                           'E_final_range': str(parameter_ranges['E']),
                                           'R_range': str(parameter_ranges['R']),
                                           'K_range': str(parameter_ranges['K']),
                                           'Mix_range': str(parameter_ranges['Mix']),
                                           'Scaling_range': str(parameter_ranges['S']),
                                           'T_range': str(parameter_ranges['T']),
                                           'n_walkers': n_walkers, 'n_steps': n_steps,
                                           'burn_in': burn_in,
                                           'csm': csm,
                                           'normalization': normalization,
                                           'Tthreshold_lum' : LumTthreshold,
                                           'time': time_now}, orient='index')
    run_param_df.to_csv(os.path.join(output_dir, 'run_parameters.csv'))



def save_3d_array(list_3d, output_dir, type, step):
    array = np.array(list_3d)
    m_shape = array.shape
    array = array.reshape(m_shape[0], m_shape[1] * m_shape[2])
    np.savetext(os.path.join(output_dir, type+ '_models_step' + str(step) + '.txt'), array)


def import_lum(SN_name):
    lum_path = os.path.join('SN_data', SN_name + '_lum')
    data_lum = pd.read_csv(lum_path)
    return data_lum


def import_veloc(SN_name):
    veloc_path = os.path.join('SN_data', SN_name + '_veloc')
    data_veloc = pd.read_csv(veloc_path)
    data_veloc.rename({'absorption_mean_velocity': 'veloc', 'absorption_std_velocity': 'dveloc'}, axis='columns',
                      inplace=True)
    if not 'dveloc' in data_veloc.columns.tolist():
        data_veloc['dveloc'] = 1.0
    data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
    data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
    return data_veloc


def import_mag(SN_name):
    mag_path = os.path.join('SN_data', SN_name + '_mag')
    data_filepath = os.path.join(mag_path)
    data_mag = pd.read_csv(data_filepath, usecols=['dmag', 'filter', 'abs_mag', 't_from_discovery'])
    filters = list(data_mag['filter'].unique())
    # filters = list(set(filters).intersection(['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']))
    filters = list(set(filters).intersection(['g', 'r', 'i', 'z', 'B', 'V', 'R', 'I']))
    data_mag = data_mag.loc[data_mag['filter'].isin(filters)]
    data_mag = data_mag.sort_values('t_from_discovery')
    return data_mag


def emcee_fit_params(SN_name, res_dir, n_walkers, n_steps, burn_in, ranges_dict, fitting_type, csm, LumTthreshold=False,
                     normalization=False, nonuniform_priors=None, init_guesses=None):
    if csm:
        n_params = 8
    else:
        n_params = 6
    if init_guesses is None:
        init_guesses = initial_guesses(ranges_dict, n_walkers, csm)
    initialize_empty_models(ranges_dict)
    data = load_SN_data(fitting_type, SN_name)
    time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    write_params_file(ranges_dict, SN_name, n_walkers, n_steps,
                                csm, LumTthreshold, normalization, burn_in, time_now, res_dir)

    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior,
                                    args=[data, ranges_dict, fitting_type, csm, LumTthreshold, normalization, nonuniform_priors])
    sampler.run_mcmc(init_guesses, n_steps)
    return sampler


def initialize_empty_models(ranges_dict):
    mod.initialize_empty_models_dict(models, ranges_dict)

