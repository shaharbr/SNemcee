import numpy as np
import emcee
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import snec_model_interpolator as interp
import pandas as pd
import os
import corner
import csv
from scipy.stats import chi2
import copy
from numpy import trapz



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
T_thresh = 10 ** 3.5
models = {}


def initialize_empty_models_dict(ranges_dict):
    for data_type in ['lum', 'veloc', 'mag', 'temp']:
        models[data_type] = {}
        for Mzams in ranges_dict['Mzams']:
            models[data_type][Mzams] = {}
            for Ni in ranges_dict['Ni']:
                models[data_type][Mzams][Ni] = {}
                for E in ranges_dict['E']:
                    models[data_type][Mzams][Ni][E] = {}
                    for R in ranges_dict['R']:
                        models[data_type][Mzams][Ni][E][R] = {}
                        for K in ranges_dict['K']:
                            models[data_type][Mzams][Ni][E][R][K] = {}
                            for Mix in ranges_dict['Mix']:
                                models[data_type][Mzams][Ni][E][R][K][Mix] = None


def get_surrouding_values(requested, ranges_dict):
    params = list(ranges_dict.keys())
    surrouding_values = {param: [] for param in params}
    for i in range(len(requested)):
        param_range = np.array(ranges_dict[params[i]])
        below = np.max(param_range[param_range <= requested[i]])
        above = np.min(param_range[param_range >= requested[i]])
        surrouding_values[params[i]] = [below, above]
    return surrouding_values


def load_model(Mzams, Ni, E, R, K, Mix, data_type):
    if R == 0 or K == 0:
        R = 0
        K = 0
    name = 'M' + str(Mzams) + \
           '_Ni' + str(Ni) + \
           '_E' + str(E) + \
           '_Mix' + str(Mix) + \
           '_R' + str(R) + \
           '_K' + str(K)
    if data_type == 'lum':
        modelpath = os.path.join('SNEC_models', name, 'lum_observed.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath)
            time_col = snec_model['t_from_discovery'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, time_col, snec_model['Lum'])
            return snec_model
    elif data_type == 'veloc':
        modelpath = os.path.join('SNEC_models', name, 'vel_Fe.dat')
        if os.stat(modelpath).st_size < 10 ** 4:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath)
            time_col = snec_model['t_from_discovery'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, time_col, snec_model['vel'])
            return snec_model
    elif data_type == 'temp':
        modelpath = os.path.join('SNEC_models', name, 'T_eff.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        snec_model = pd.read_csv(modelpath)
        time_col = snec_model['t_from_discovery'] / 86400  # sec to days
        interp_days = np.linspace(0, 200, 2001)
        snec_model = np.interp(interp_days, time_col, snec_model['temp'])
        return snec_model
    elif data_type == 'mag':
        modelpath = os.path.join('SNEC_models', name, 'magnitudes_pys.dat')
        lumpath = os.path.join('SNEC_models', name, 'lum_observed.dat')
        if os.stat(lumpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            mag_file = pd.read_csv(modelpath, names=['t_from_discovery', 'Teff', 'PTF_R_AB', 'R^2', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'])
            time_col = mag_file['t_from_discovery'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model_dict = {}
            for filter in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
                snec_model_dict[filter] = np.interp(interp_days, time_col, mag_file[filter])
            snec_model_dict['t_from_discovery'] = interp_days
            snec_model = pd.DataFrame(snec_model_dict)
            snec_model = snec_model.sort_values('t_from_discovery')
            return snec_model


def load_surrounding_models(requested, ranges_dict, fitting_type, Tthreshold):
    any_Tthreshold = any([Tthreshold['lum'], Tthreshold['veloc'], Tthreshold['mag']])
    surrouding_values = get_surrouding_values(requested, ranges_dict)
    for Mzams in surrouding_values['Mzams']:
        for Ni in surrouding_values['Ni']:
            for E in surrouding_values['E']:
                for R in surrouding_values['R']:
                    for K in surrouding_values['K']:
                        for Mix in surrouding_values['Mix']:
                            if 'lum' in fitting_type or 'combined' in fitting_type:
                                if models['lum'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['lum'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'lum')
                            if 'veloc' in fitting_type or 'combined' in fitting_type:
                                if models['veloc'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['veloc'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'veloc')
                            if 'mag' in fitting_type or 'combined' in fitting_type:
                                if models['mag'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['mag'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'mag')
                            if any_Tthreshold:
                                if models['temp'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['temp'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'temp')


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

def chi_square_norm(y, dy, y_fit):
    return np.sum(((y - y_fit) / dy) ** 2)


# TODO this step can introduce some bias - SNe that don't have early velocities
#  will select against models that cool fast
def temp_thresh_cutoff(requested_theta, surrounding_values, models, data_x):
    temp_fit = interp.snec_interpolator(requested_theta[0:6], surrounding_values, models['temp'],
                                        data_x)
    max_x = 0
    if not isinstance(temp_fit, str):
        if len(temp_fit[temp_fit >= T_thresh]) > 0:
            min_temp_above_Tthresh = np.min(temp_fit[temp_fit >= T_thresh])
            x_out = data_x[temp_fit > min_temp_above_Tthresh]
            max_x = np.max(x_out)
    return max_x

def calc_lum_likelihood(theta, data_dict, surrounding_values, Tthreshold_dict, normalization):
    data = data_dict['lum']
    Tthreshold = Tthreshold_dict['lum']
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    if Tthreshold:
        max_x = temp_thresh_cutoff(theta[0:6], surrounding_values, models, data_x_moved)
        data_x_moved = data_x_moved[data_x_moved <= max_x]
        data = data.loc[data['t_from_discovery'] - theta[7] <= max_x]
    data_y = data['Lum']
    data_dy = data['dLum0']
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['lum'], data_x_moved)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        # calculate the log likelihood
        df = len(data_y) - 1
        if normalization:
            norm = len(data_y)
        else:
            norm = 1
        return chi2.logpdf(chi_square_norm(data_y, data_dy, y_fit), df) / norm
    else:
        # print('impossible SN')
        return - np.inf


def calc_veloc_likelihood(theta, data_dict, surrounding_values, Tthreshold_dict, normalization):
    data = data_dict['veloc']
    Tthreshold = Tthreshold_dict['veloc']
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    if Tthreshold:
        max_x = temp_thresh_cutoff(theta[0:6], surrounding_values, models, data_x_moved)
        data_x_moved = data_x_moved[data_x_moved <= max_x]
        data = data.loc[data['t_from_discovery'] - theta[7] <= max_x]
    data_y = data['veloc']
    data_dy = data['dveloc']
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['veloc'], data_x_moved)
    if not isinstance(y_fit, str):
        # calculate the log likelihood
        df = len(data_y) - 1
        if normalization:
            norm = len(data_y)
        else:
            norm = 1
        return chi2.logpdf(chi_square_norm(data_y, data_dy, y_fit), df) / norm
    else:
        # print('impossible SN')
        return - np.inf


def calc_mag_likelihood(theta, data_dict, surrounding_values, Tthreshold_dict, normalization):
    data = data_dict['mag']
    Tthreshold = Tthreshold_dict['mag']
    log_likeli = 0
    any_filter_data = False
    y_fit = {}
    filters = list(data['filter'].unique())
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    if Tthreshold:
        max_x = temp_thresh_cutoff(theta[0:6], surrounding_values, models, data_x_moved)
        data = data.loc[data['t_from_discovery'] - theta[7] <= max_x]
    for filt in filters:
        data_filt = data.loc[data['filter'] == filt]
        data_x_filt_moved = data_filt['t_from_discovery'] - theta[7]
        data_y_filt = data_filt['abs_mag']
        data_dy_filt = data_filt['dmag']
        y_fit[filt] = interp.snec_interpolator(theta[0:6], surrounding_values, models['mag'], data_x_filt_moved,
                                               filter=filt)
        if not isinstance(y_fit[filt], str):
            # multiply whole graph by scaling factor
            y_fit[filt] = y_fit[filt] * theta[6]
            # calculate the log likelihood
            df = len(data_y_filt) - 1
            if normalization:
                norm = len(data_y_filt)
            else:
                norm = 1
            log_likeli += chi2.logpdf(chi_square_norm(data_y_filt, data_dy_filt, y_fit[filt]), df) / norm
            any_filter_data = True
    if any_filter_data:
        return log_likeli
    else:
        return - np.inf


def log_likelihood(theta, data, ranges_dict, fitting_type, Tthreshold, normalization):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form described above

    fitting_type : lum, mag, veloc, combined or combined_normalized
    """
    # print(theta)
    # ranges_list = dict_to_list(ranges_dict)
    # print(ranges_list)
    if theta_in_range(theta, ranges_dict):
        # print('ok SN')
        surrounding_values = get_surrouding_values(theta[0:6], ranges_dict)
        load_surrounding_models(theta[0:6], ranges_dict, fitting_type, Tthreshold)
        args = [theta, data, surrounding_values, Tthreshold, normalization]
        if fitting_type == 'lum':
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


def log_posterior(theta, data, ranges_dict, fitting_type, csm, Tthreshold, normalization, nonuniform_priors):
    if not csm:
        theta =list(np.insert(theta, 3, np.zeros((2,)), axis=0))
    if theta_in_range(theta, ranges_dict):
        lp = log_prior(theta, ranges_dict, nonuniform_priors)
        ll = log_likelihood(theta, data, ranges_dict, fitting_type, Tthreshold, normalization)
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


def emcee_fit_params(data, n_walkers, n_steps, ranges_dict, fitting_type, csm, Tthreshold, normalization=False, nonuniform_priors=None, init_guesses=None):
    if csm:
        n_params = 8
    else:
        n_params = 6
    if init_guesses is None:
        init_guesses = initial_guesses(ranges_dict, n_walkers, csm)
    initialize_empty_models_dict(ranges_dict)
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior,
                                    args=[data, ranges_dict, fitting_type, csm, Tthreshold, normalization, nonuniform_priors])
    sampler.run_mcmc(init_guesses, n_steps)
    return sampler


def corner_plot(flat_sampler_no_burnin, ranges_dict, output_dir):
    labels = list(ranges_dict.keys())
    corner_range = [1.] * len(labels)
    f_corner = corner.corner(flat_sampler_no_burnin, labels=labels, range=corner_range)
    f_corner.savefig(os.path.join(output_dir, 'corner_plot.png'))


def chain_plots(sampler_chain, ranges_dict, output_dir, burn_in, first_stage_steps=None):
    keys = list(ranges_dict.keys())
    for i in range(len(keys)):
        key = keys[i]
        plt.figure()
        plt.plot(sampler_chain[:, :, i].T, color='k', alpha=0.1)
        plt.xlabel('Step Number')
        plt.ylabel(key)
        plt.axvspan(0, burn_in, alpha=0.1, color='grey')
        if first_stage_steps is not None:
            plt.axvline(x=first_stage_steps, color='black')
            plt.axvspan(first_stage_steps, first_stage_steps+burn_in, alpha=0.1, color='grey')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, key+'.png'))


def get_each_walker_result(sampler_chain, ranges_dict, step):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
        last_results = sampler_chain[:, step:, i]
        avg = np.mean(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg


def get_param_results_dict(sampler_chain, ranges_dict, step, output_dir):
    params = list(ranges_dict.keys())
    dict = {}
    for i, param in enumerate(params):
        last_results = sampler_chain[:, step:, i]
        avg = np.mean(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[param] = avg
        dict[param +'_lower'] = avg - sigma_lower
        dict[param + '_upper'] = sigma_upper - avg
    with open(os.path.join(output_dir, 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)
    return dict


def write_params_file(parameter_ranges, S_range, SN_name, n_walkers, n_steps, csm, Tthreshold_dict, normalization, burn_in, time_now, output_dir):
    run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                                           'Mzams_range': str(parameter_ranges['Mzams']),
                                           'Ni_range': str(parameter_ranges['Ni']),
                                           'E_final_range': str(parameter_ranges['E']),
                                           'R_range': str(parameter_ranges['R']),
                                           'K_range': str(parameter_ranges['K']),
                                           'Mix_range': str(parameter_ranges['Mix']),
                                           'Scaling_range': str(S_range),
                                           'T_range': str(parameter_ranges['T']),
                                           'n_walkers': n_walkers, 'n_steps': n_steps,
                                           'burn_in': burn_in,
                                           'csm': csm,
                                           'normalization': normalization,
                                           'Tthreshold_lum' : Tthreshold_dict['lum'],
                                           'Tthreshold_mag': Tthreshold_dict['mag'],
                                           'Tthreshold_veloc': Tthreshold_dict['veloc'],
                                           'time': time_now}, orient='index')
    run_param_df.to_csv(os.path.join(output_dir, 'run_parameters.csv'))


def rounded_str(x):
    if not np.isinf(x) and np.abs(x) > 0.0000001:
        rounded = round(x, 2-int(np.floor(np.log10(abs(x)))))
        if rounded > 100:
            rounded = int(rounded)
    else:
        rounded = 0
    return str(rounded)


def result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir):
    param_dict = get_param_results_dict(sampler_chain, ranges_dict, step, output_dir)
    res_text = SN_name + '\n'
    params = list(ranges_dict.keys())
    for param in params:
        if (param != 'K') & (param != 'R'):
            res_text += param + ': ' + rounded_str(param_dict[param]) + r'$\pm$ [' +\
                            rounded_str(param_dict[param+'_lower']) + ',' +\
                            rounded_str(param_dict[param+'_upper']) + ']\n'
    if (param_dict['K'] == 0) & (param_dict['R'] == 0):
        res_text += 'no CSM'
    else:
        res_text += 'K' + ': ' + rounded_str(param_dict['K']) + r'$\pm$ [' + \
                    rounded_str(param_dict['K_lower']) + ',' + \
                    rounded_str(param_dict['K_upper']) + ']\n'
        res_text += 'R' + ': ' + rounded_str(param_dict['R']) + r'$\pm$ [' + \
                    rounded_str(param_dict['R_lower']) + ',' + \
                    rounded_str(param_dict['R_upper']) + ']\n'
    return res_text


def plot_lum_with_fit(data_dict, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, Tthreshold_dict, normalization):
    data = data_dict['lum']
    Tthreshold = Tthreshold_dict['lum']
    data_x = data['t_from_discovery']
    data_y = data['Lum']
    dy0 = data['dLum0']
    dy1 = data['dLum1']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict):
            load_surrounding_models(requested[0:6], ranges_dict, 'lum', Tthreshold_dict)
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            data_x_moved = data_x - T
            if Tthreshold:
                max_x = temp_thresh_cutoff(requested, surrounding_values, models, data_x_moved)
                x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            else:
                x_plotting = np.linspace(-T, 200-T, 2001)
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'], x_plotting)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                ax.plot(x_plotting, np.log10(y_fit), alpha=0.4)
                log_likeli.append(calc_lum_likelihood(requested, data_dict, surrounding_values, Tthreshold_dict, normalization))
    log_likeli = np.mean(log_likeli)
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 200)
    ax.set_ylim(39.5, 43.5)
    ax.set_title('step ' + str(step), fontsize=14)
                 # + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + '.png'))
    return log_likeli


def plot_mag_with_fit(data_dict, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, Tthreshold_dict, normalization):
    data = data_dict['mag']
    Tthreshold = Tthreshold_dict['mag']
    filters = list(data['filter'].unique())
    data_x = data['t_from_discovery']
    y_fit = {}
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict):
            load_surrounding_models(requested[0:6], ranges_dict, 'mag', Tthreshold_dict)
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            data_x_moved = data_x - T
            if Tthreshold:
                max_x = temp_thresh_cutoff(requested, surrounding_values, models, data_x_moved)
                x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            else:
                x_plotting = np.linspace(-T, 200-T, 2001)
            for filt in filters:
                y_fit[filt] = interp.snec_interpolator(requested[0:6], surrounding_values,
                                                       models['mag'], x_plotting, filter=filt)
                if not isinstance(y_fit, str):
                    # multiply whole graph by scaling factor
                    y_fit[filt] = y_fit[filt] * S
                    ax.plot(x_plotting, y_fit[filt], color=colors[filt], alpha=0.3)
                    log_likeli.append(calc_mag_likelihood(requested, data_dict, surrounding_values, Tthreshold_dict, normalization))
    log_likeli = np.mean(log_likeli)
    for filt in filters:
        data_filt = data.loc[data['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_y = data_filt['abs_mag']
        data_dy = data_filt['dmag']
        ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None',
                    label=filt, color=colors[filt])
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_title('step ' + str(step), fontsize=14)
                 # + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    # TODO there was an error here that sometimes didn't allow printing the log likelihood title
    plt.tight_layout()
    ax.set_xlim(-2, 137)
    ax.invert_yaxis()
    ax.set_ylim(-14, -19)
    f_fit.savefig(os.path.join(output_dir, 'mag_lightcurve_fit' + str(step) + '.png'))
    return log_likeli


def plot_veloc_with_fit(data_dict, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, Tthreshold_dict, normalization):
    data = data_dict['veloc']
    Tthreshold = Tthreshold_dict['veloc']
    data_x = data['t_from_discovery']
    data_y = data['veloc']
    data_dy = data['dveloc']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict):
            load_surrounding_models(requested[0:6], ranges_dict, 'veloc', Tthreshold_dict)
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            data_x_moved = data_x - T
            if Tthreshold:
                max_x = temp_thresh_cutoff(requested, surrounding_values, models, data_x_moved)
                x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            else:
                x_plotting = np.linspace(-T, 150 - T, 1501)
            x_plotting_moved = x_plotting - T
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['veloc'], x_plotting_moved)
            if not isinstance(y_fit, str):
                ax.plot(x_plotting, y_fit, alpha=0.4)
                log_likeli.append(calc_veloc_likelihood(requested, data_dict, surrounding_values, Tthreshold_dict, normalization))
    log_likeli = np.mean(log_likeli)
    ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 150)
    ax.set_title('step '+str(step), fontsize=14)
                 # +'\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, 'velocities_fit' +str(step) +'.png'))
    return log_likeli


def plot_fit_with_data(sampler_chain, data, ranges_dict, fitting_type, output_dir, SN_name, step, Tthreshold, normalization=False):
    n_walkers, n_steps, n_params = np.shape(sampler_chain)
    args = (data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, Tthreshold, normalization)
    if fitting_type == 'lum':
        log_likeli = plot_lum_with_fit(*args)
    if fitting_type == 'mag':
        log_likeli = plot_mag_with_fit(*args)
    if fitting_type == 'veloc':
        log_likeli = plot_veloc_with_fit(*args)
    if fitting_type == 'lum-veloc':
        log_likeli = plot_lum_with_fit(*args) + \
                     plot_veloc_with_fit(*args)
    if fitting_type == 'mag-veloc':
        log_likeli = plot_mag_with_fit(*args) + \
                     plot_veloc_with_fit(*args)
    if fitting_type == 'combined':
        log_likeli = plot_lum_with_fit(*args) + \
                     plot_mag_with_fit(*args)+ \
                     plot_veloc_with_fit(*args)
    return log_likeli


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
    filters = list(set(filters).intersection(['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']))
    data_mag = data_mag.loc[data_mag['filter'].isin(filters)]
    data_mag = data_mag.sort_values('t_from_discovery')
    return data_mag