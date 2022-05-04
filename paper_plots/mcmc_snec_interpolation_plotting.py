import emcee
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import snec_model_interpolator_plotting as interp
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
T_thresh = 10 ** 3.75
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


def load_model(Mzams, Ni, E, R, K, Mix, data_type, extend_tail=False):
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
        modelpath = os.path.join('..', 'all_lum_data', name, 'lum_observed.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath,
                                     names=['t_from_discovery', 'Lum'], sep=r'\s+')
            time_col = snec_model['t_from_discovery'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, time_col, snec_model['Lum'])
            if extend_tail is not False:
                last_30d_x = interp_days[-100:]
                last_30d_y = snec_model[-100:]
                last_30d_ylog = np.log(last_30d_y)
                tail_poly1d = np.poly1d(np.polyfit(last_30d_x, last_30d_ylog, deg=1))
                extension_days = np.linspace(200.1, 200+extend_tail, int(10*extend_tail))
                extension_lumlog = np.array([tail_poly1d(extension_days[i]) for i in range(len(extension_days))])
                extension_lum = np.exp(extension_lumlog)
                snec_model = np.concatenate((snec_model, extension_lum))
            return snec_model
    elif data_type == 'veloc':
        modelpath = os.path.join('..', 'all_veloc_data', name, 'vel_Fe.dat')
        if os.stat(modelpath).st_size < 10 ** 4:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath)
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, snec_model['t_from_discovery'],
                                   snec_model['veloc'])
            return snec_model
    elif data_type == 'temp':
        modelpath = os.path.join('..', 'all_temp_rad_data', name, 'T_eff.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        snec_model = pd.read_csv(modelpath,
                                 names=['t_from_discovery', 'temp'], sep=r'\s+')
        time_col = snec_model['t_from_discovery'] / 86400  # sec to days
        interp_days = np.linspace(0, 200, 2001)
        snec_model = np.interp(interp_days, time_col, snec_model['temp'])
        return snec_model
    elif data_type == 'mag':
        modelpath = os.path.join('..', 'all_pys_mag_data', name, 'magnitudes.dat')
        lumpath = os.path.join('..', 'all_lum_data', name, 'lum_observed.dat')
        if os.stat(lumpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            mag_file = pd.read_csv(modelpath,
                                   names=['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'],
                                   sep=r'\s+')
            mag_file = mag_file.abs()
            time_col = mag_file['time'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model_dict = {}
            for filter in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
                snec_model_dict[filter] = np.interp(interp_days, time_col, mag_file[filter])
            snec_model_dict['time'] = interp_days
            snec_model = pd.DataFrame(snec_model_dict)
            snec_model = snec_model.sort_values('time')
            return snec_model


def load_surrounding_models(requested, ranges_dict, fitting_type, extend_tail=False):
    surrouding_values = get_surrouding_values(requested, ranges_dict)
    for Mzams in surrouding_values['Mzams']:
        for Ni in surrouding_values['Ni']:
            for E in surrouding_values['E']:
                for R in surrouding_values['R']:
                    for K in surrouding_values['K']:
                        for Mix in surrouding_values['Mix']:
                            if 'lum' in fitting_type:
                                if models['lum'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['lum'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'lum', extend_tail)
                            if 'veloc' in fitting_type:
                                if models['veloc'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['veloc'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'veloc')
                                if models['temp'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['temp'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'temp')
                            if 'mag' in fitting_type:
                                if models['mag'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['mag'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'mag')


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
                    print('logprior', log_prior)
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


def calc_lum_likelihood(theta, data, surrounding_values, extend_tail=False):
    data_x_moved = data['t_from_discovery'] - theta[7]
    data_y = data['Lum']
    data_dy = data['dLum0']
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['lum'], data_x_moved, extend_tail)
    exit()
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        # calculate the log likelihood
        df = len(data_y) - 1
        print(chi2.logpdf(chi_square_norm(data_y, data_dy, y_fit), df))
        return chi2.logpdf(chi_square_norm(data_y, data_dy, y_fit), df)
    else:
        print('impossible SN')
        return - np.inf

def calc_veloc_likelihood(theta, data, surrounding_values, Tthreshold=False):
    data_y = data['veloc']
    data_dy = data['dveloc']
    data_x = data['t_from_discovery']
    data_x_moved = data_x - theta[7]
    if Tthreshold:
        temp_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['temp'], data_x_moved)
        if not isinstance(temp_fit, str):
            below_Tthresh = temp_fit[temp_fit <= T_thresh]
            print('tthreshhhhhhhhhhh')
            print(below_Tthresh)
            print(len(below_Tthresh))
            if len(below_Tthresh) > 0:
                max_temp_below_Tthresh = np.max(below_Tthresh)
                data_x_moved = data_x_moved[temp_fit > max_temp_below_Tthresh]
                if len(data_x_moved) <= 1:
                    # TODO this step introduces some bias - SNe that don't have early velocities
                    #  will select against models that cool fast
                    print('cooled too fast, no early velocity data')
                    return - np.inf
                data = data.loc[temp_fit > max_temp_below_Tthresh]
                data_y = data['veloc']
                data_dy = data['dveloc']
            else:
                print('cooled too fast')
                return - np.inf
        else:
            print('impossible SN')
            return - np.inf
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['veloc'], data_x_moved)
    if not isinstance(y_fit, str):
        # calculate the log likelihood
        df = len(data_y) - 1
        return chi2.logpdf(chi_square_norm(data_y, data_dy, y_fit), df)
    else:
        print('impossible SN')
        return - np.inf


def calc_mag_likelihood(theta, data, surrounding_values):
    log_likeli = 0
    data_x_allfilt = data['t_from_discovery'] - theta[7]
    filters = list(data['filter'].unique())
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['mag'], data_x_allfilt)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        y_fit['time'] = data_x_allfilt
        for filt in filters:
            data_filt = data.loc[data['filter'] == filt]
            data_x_filt = data_filt['t_from_discovery'] - theta[7]
            data_y_filt = data_filt['abs_mag']
            data_dy_filt = data_filt['dmag']
            y_fit_filt = y_fit[filt].loc[y_fit['time'].isin(data_x_filt)]
            # calculate the log likelihood
            df = len(data_y_filt) - 1
            log_likeli += chi2.logpdf(chi_square_norm(data_y_filt, data_dy_filt, y_fit_filt), df)
            # print('log likelihood', log_likeli)
        return log_likeli
    else:
        print('impossible SN')
        return - np.inf


def log_likelihood(theta, data, ranges_dict, fitting_type, extend_tail=False):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form described above

    fitting_type : lum, mag, veloc, combined or combined_normalized
    """
    print(theta)
    # ranges_list = dict_to_list(ranges_dict)
    # print(ranges_list)
    if theta_in_range(theta, ranges_dict) or theta[3] == 0:
        print('ok SN')
        surrounding_values = get_surrouding_values(theta[0:6], ranges_dict)
        load_surrounding_models(theta[0:6], ranges_dict, fitting_type, extend_tail)
        if fitting_type == 'lum':
            log_likeli = calc_lum_likelihood(theta, data['lum'], surrounding_values, extend_tail)
        elif fitting_type == 'mag':
            log_likeli = calc_mag_likelihood(theta, data['mag'], surrounding_values)
        elif fitting_type == 'veloc':
            log_likeli = calc_veloc_likelihood(theta, data['veloc'], surrounding_values)
        elif fitting_type == 'lum_veloc':
            log_likeli = calc_lum_likelihood(theta, data['lum'], surrounding_values, extend_tail) + \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values)
        elif fitting_type == 'lum_veloc_Tthresh_normalized':
            log_likeli = calc_lum_likelihood(theta, data['lum'], surrounding_values, extend_tail)+ \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values, Tthreshold=True)
        elif fitting_type == 'lum_veloc_normalized':
            num_obs_lum = len(data['lum']['t_from_discovery'])
            num_obs_veloc = len(data['veloc']['t_from_discovery'])
            log_likeli = calc_lum_likelihood(theta, data['lum'], surrounding_values, extend_tail) / num_obs_lum + \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values) / num_obs_veloc
        elif fitting_type == 'mag_veloc':
            log_likeli = calc_mag_likelihood(theta, data['mag'], surrounding_values) + \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values)
        elif fitting_type == 'mag_veloc_normalized':
            num_obs_mag = len(data['mag']['t_from_discovery'])
            num_obs_veloc = len(data['veloc']['t_from_discovery'])
            log_likeli = calc_mag_likelihood(theta, data['mag'], surrounding_values) / num_obs_mag + \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values) / num_obs_veloc
        elif fitting_type == 'combined':
            log_likeli = calc_lum_likelihood(theta, data['lum'], surrounding_values, extend_tail) + \
                         calc_mag_likelihood(theta, data['mag'], surrounding_values) + \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values)
        elif fitting_type == 'combined_normalized':
            num_obs_lum = len(data['lum']['t_from_discovery'])
            num_obs_mag = len(data['mag']['t_from_discovery'])
            num_obs_veloc = len(data['veloc']['t_from_discovery'])
            log_likeli = calc_lum_likelihood(theta, data['lum'], surrounding_values, extend_tail) / num_obs_lum+ \
                         calc_mag_likelihood(theta, data['mag'], surrounding_values) / num_obs_mag+ \
                         calc_veloc_likelihood(theta, data['veloc'], surrounding_values) / num_obs_veloc
        else:
            print('fitting_type should be: lum, mag, veloc, lum_veloc, lum_veloc_normalized, mag_veloc, mag_veloc_normalized, combined or combined_normalized')
    else:
        print('log lik')
        print(theta)
        print(ranges_dict)
        print('out of range')
        return - np.inf  # just a very big number so it won't go past the edge values
    print('loglik', log_likeli)
    return log_likeli


def log_posterior(theta, data, ranges_dict, fitting_type, nonuniform_priors=None, extend_tail=False):
    if theta_in_range(theta, ranges_dict):
        lp = log_prior(theta, ranges_dict, nonuniform_priors)
        ll = log_likelihood(theta, data, ranges_dict, fitting_type, extend_tail)
        log_post = lp + ll
        print('logpost', log_post)
        return log_post
    else:
        print('log pos')
        print(theta)
        print(ranges_dict)
        print('out of range')
        return - np.inf  # just a very big number so it won't go past the edge values


def initial_guesses(ranges_dict, n_walkers):
    guesses = []
    for param in list(ranges_dict.keys()):
        guesses.append(np.random.rand(n_walkers) *
                       (ranges_dict[param][-1] - ranges_dict[param][0]) + ranges_dict[param][0])
    return np.array(guesses).T


def emcee_fit_params(data, n_walkers, n_steps, ranges_dict, fitting_type, nonuniform_priors=None, init_guesses=None, extend_tail=False):
    data_lum = data['lum']
    data_no_early_lum = copy.deepcopy(data)
    data_no_early_lum['lum'] = data_lum.loc[data_lum['t_from_discovery'] > 30]
    n_params = len(ranges_dict.keys())
    if init_guesses is None:
        init_guesses = initial_guesses(ranges_dict, n_walkers)
    initialize_empty_models_dict(ranges_dict)
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior,
                                    args=[data_no_early_lum, ranges_dict, fitting_type, nonuniform_priors, extend_tail])
    sampler.run_mcmc(init_guesses, n_steps)
    return sampler


def corner_plot(sampler_chain_flat, ranges_dict, output_dir):
    labels = list(ranges_dict.keys())
    corner_range = [1.] * len(labels)
    f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
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
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, key+'.png'))


def get_each_walker_result(sampler_chain, ranges_dict, step):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
        last_results = sampler_chain[:, step:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg


def get_param_results_dict(sampler_chain, ranges_dict, step, output_dir):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
        last_results = sampler_chain[:, step:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        # dict[params[i]+'_all_walkers'] = last_results
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    with open(os.path.join(output_dir, 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)
    return dict

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

# TODO fix the interp here to the updated

def plot_lum_with_fit(data_lum, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail=False):
    data_x = data_lum['t_from_discovery']
    data_y = data_lum['Lum']
    dy0 = data_lum['dLum0']
    dy1 = data_lum['dLum1']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    if extend_tail is not False:
        x_plotting = np.linspace(0, 200+extend_tail, int(1+10*200+extend_tail))
    else:
        x_plotting = np.linspace(0, 200, 2001)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved = x_plotting - T
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            load_surrounding_models(requested[0:6], ranges_dict, 'lum', extend_tail)
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'], data_x_moved, extend_tail)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                ax.plot(x_plotting, np.log10(y_fit), alpha=0.4)
                log_likeli.append(calc_lum_likelihood(requested, data_lum, surrounding_values, extend_tail))
    log_likeli = np.average(log_likeli)
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    if extend_tail is not False:
        ax.set_xlim(-2, 200+extend_tail)
    else:
        ax.set_xlim(-2, 200)
    ax.set_ylim(40.8, 43)
    # ax.set_ylim(np.min(data_y) * 0.5, np.max(data_y) * 5)
    ax.set_title('step ' + str(step)
                 + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    # ax.set_ylim(float(3.0 * 10 ** 39), float(1.6 * 10 ** 43))
    # ax.set_yscale('log')
    # f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + 'log.png'))
    # ax.set_yscale('linear')
    f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + '.png'))
    return log_likeli

def plot_mag_with_fit(data_mag, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name):
    filters = list(data_mag['filter'].unique())
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    x_plotting = np.linspace(0, 200, 2001)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict) or R == 0:
            load_surrounding_models(requested[0:6], ranges_dict, 'mag')
            data_x_moved_all = x_plotting - T
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['mag'], data_x_moved_all)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                for filt in filters:
                    ax.plot(x_plotting, y_fit[filt], color=colors[filt], alpha=0.3)
                log_likeli.append(calc_mag_likelihood(requested, data_mag, surrounding_values))
    log_likeli = np.average(log_likeli)
    for filt in filters:
        data_filt = data_mag.loc[data_mag['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_y = data_filt['abs_mag']
        data_dy = data_filt['dmag']
        ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None',
                    label=filt, color=colors[filt])
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_title('step ' + str(step)
                 + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    ax.set_xlim(-2, 137)
    ax.set_ylim(14, 20)
    f_fit.savefig(os.path.join(output_dir, 'mag_lightcurve_fit' + str(step) + '.png'))
    return log_likeli



def plot_veloc_with_fit(data_veloc, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, Tthreshold=False):
    data_x = data_veloc['t_from_discovery']
    data_y = data_veloc['veloc']
    data_dy = data_veloc['dveloc']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    x_plotting = np.linspace(0, np.max(data_veloc['t_from_discovery']), int(1+10*np.max(data_veloc['t_from_discovery'])))
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict) or R == 0:
            load_surrounding_models(requested[0:6], ranges_dict, 'veloc')
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            temp_x_moved = data_x - T
            if Tthreshold:
                temp_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['temp'], temp_x_moved)
                if not isinstance(temp_fit, str):
                    max_temp_below_Tthresh = np.max(temp_fit[temp_fit <= T_thresh])
                    temp_x_moved = temp_x_moved[temp_fit > max_temp_below_Tthresh]
                    if len(temp_x_moved) <= 1:
                        print('cooled too fast, no early velocity data')
                        x_plotting = []
                    else:
                        max_x_moved = np.max(temp_x_moved)
                        x_plotting = np.linspace(0, max_x_moved, int(1+max_x_moved*10))
            else:
                x_plotting = np.linspace(0, 200, 2001)
            if len(x_plotting) > 0:
                data_x_moved = x_plotting - T
                y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['veloc'], data_x_moved)
                if not isinstance(y_fit, str):
                    ax.plot(x_plotting, y_fit, alpha=0.4)
                    log_likeli.append(calc_veloc_likelihood(requested, data_veloc, surrounding_values, Tthreshold))
    log_likeli = np.average(log_likeli)
    ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 140)
    ax.set_title('step '+str(step)
                 +'\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, 'velocities_fit' +str(step) +'.png'))
    return log_likeli



def plot_lightcurve_with_fit(sampler_chain, SN_data, ranges_dict, fitting_type, output_dir, n_walkers, SN_name, step, extend_tail=False):
    log_likeli = 0
    if fitting_type == 'lum':
        log_likeli = plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail)
    if fitting_type == 'mag':
        log_likeli = plot_mag_with_fit(SN_data['mag'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'veloc':
        log_likeli = plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'lum_veloc':
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail)
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'lum_veloc_Tthresh_normalized':
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail)
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, Tthreshold=True)
    if fitting_type == 'lum_veloc_normalized':
        log_likeli = 0
        num_obs_veloc = len(SN_data['veloc']['t_from_discovery'])
        num_obs_lum = len(SN_data['lum']['t_from_discovery'])
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail) / num_obs_lum
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_veloc
    if fitting_type == 'mag_veloc':
        log_likeli = 0
        log_likeli += plot_mag_with_fit(SN_data['mag'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'mag_veloc_normalized':
        log_likeli = 0
        num_obs_veloc = len(SN_data['veloc']['t_from_discovery'])
        num_obs_mag = len(SN_data['mag']['t_from_discovery'])
        log_likeli += plot_mag_with_fit(SN_data['mag'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_mag
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_veloc
    if fitting_type == 'combined':
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail)
        log_likeli += plot_mag_with_fit(SN_data['mag'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'combined_normalized':
        log_likeli = 0
        num_obs_mag = len(SN_data['mag']['t_from_discovery'])
        num_obs_veloc = len(SN_data['veloc']['t_from_discovery'])
        num_obs_lum = len(SN_data['lum']['t_from_discovery'])
        log_likeli += plot_mag_with_fit(SN_data['mag'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_mag
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, extend_tail) / num_obs_lum
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_veloc
    return log_likeli


