import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import os
import mcmc_snec
import re
import snec_model_interpolator as interp
import corner
from scipy.stats import chi2
import colorful_corner_plot as color_corner


T_thresh = 10 ** 3.5
extend_tail = False
models = {}

colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}

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
    for i, req in enumerate(requested):
        param_range = np.array(ranges_dict[params[i]])
        below = np.max(param_range[param_range <= req])
        above = np.min(param_range[param_range >= req])
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


def import_ranges(list, run_params):
    ranges = []
    for name in list:
        if 'range' in name:
            content = re.sub(r'[\[\]\s]', '', run_params.iloc[0][name]).split(',')
            if 'R_' in name or 'K_' in name:
                content = [int(content[i]) for i in range(len(content))]
            else:
                content = [float(content[i]) for i in range(len(content))]
            ranges.append(content)
    ranges_dict = {'Mzams': ranges[0], 'Ni': ranges[1], 'E': ranges[2],
                   'R': ranges[3], 'K': ranges[4], 'Mix': ranges[5],
                   'S': ranges[6], 'T': ranges[7]}
    return ranges_dict


def get_param_results_dict(sampler_df, ranges_dict):
    params = list(ranges_dict.keys())
    dict = {}
    for param in params:
        avg = np.average(sampler_df[param])
        sigma_lower, sigma_upper = np.percentile(sampler_df[param], [16, 84])
        dict[param] = avg
        dict[param + '_lower'] = avg - sigma_lower
        dict[param + '_upper'] = sigma_upper - avg
    return dict


def result_text_from_dict(sampler_df, ranges_dict):
    param_dict = get_param_results_dict(sampler_df, ranges_dict)
    params = list(ranges_dict.keys())
    res_text = ''
    for param in params:
        if (param != 'K') & (param != 'R'):
            res_text += param + ': ' + mcmc_snec.rounded_str(param_dict[param]) + r'$\pm$ [' +\
                            mcmc_snec.rounded_str(param_dict[param+'_lower']) + ',' +\
                            mcmc_snec.rounded_str(param_dict[param+'_upper']) + ']\n'
    if (param_dict['K'] == 0) & (param_dict['R'] == 0):
        res_text += 'no CSM'
    else:
        res_text += 'K' + ': ' + mcmc_snec.rounded_str(param_dict['K']) + r'$\pm$ [' + \
                    mcmc_snec.rounded_str(param_dict['K_lower']) + ',' + \
                    mcmc_snec.rounded_str(param_dict['K_upper']) + ']\n'
        res_text += 'R' + ': ' + mcmc_snec.rounded_str(param_dict['R']) + r'$\pm$ [' + \
                    mcmc_snec.rounded_str(param_dict['R_lower']) + ',' + \
                    mcmc_snec.rounded_str(param_dict['R_upper']) + ']\n'
    return res_text


# def chi_square_norm(y, dy, y_fit, df):
#     return np.sum(((y - y_fit) / dy) ** 2) / df


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


def calc_lum_likelihood(theta, data, surrounding_values, Tthreshold_dict, normalization):
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
        return chi2.logpdf(mcmc_snec.chi_square_norm(data_y, data_dy, y_fit), df) / norm
    else:
        print('impossible SN')
        return - np.inf


def calc_veloc_likelihood(theta, data, surrounding_values, Tthreshold_dict, normalization):
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
        return chi2.logpdf(mcmc_snec.chi_square_norm(data_y, data_dy, y_fit), df) / norm
    else:
        print('impossible SN')
        return - np.inf


def calc_mag_likelihood(theta, data, surrounding_values, Tthreshold_dict, normalization):
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
            log_likeli += chi2.logpdf(mcmc_snec.chi_square_norm(data_y_filt, data_dy_filt, y_fit[filt]), df) / norm
            any_filter_data = True
    if any_filter_data:
        return log_likeli
    else:
        return - np.inf


def plot_lum_with_fit(data, sampler_df, ranges_dict, n_walkers, ax, Tthreshold_dict, normalization):
    print('Tthreshold_dict')
    print(Tthreshold_dict)
    Tthreshold = Tthreshold_dict['lum']
    data_x = data['t_from_discovery']
    data_y = data['Lum']
    dy0 = data['dLum0']
    dy1 = data['dLum1']
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            load_surrounding_models(requested[0:6], ranges_dict, 'lum', Tthreshold_dict)
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            data_x_moved = data_x - T
            if Tthreshold:
                max_x = temp_thresh_cutoff(requested, surrounding_values, models, data_x_moved)
                print('Tthreshold activ', max_x)
                x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            else:
                x_plotting = np.linspace(-T, 200-T, 2001)
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'], x_plotting)
            # y_fit_on_data_times = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'],
            #                                                data_x_moved)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                # y_fit_on_data_times = y_fit_on_data_times * S
                ax.plot(x_plotting, np.log10(y_fit), alpha=0.1, color='purple')
                print('lum')
                print(data_x)
                print(max_x)
                print(x_plotting)
                log_likeli.append(calc_lum_likelihood(requested, data, surrounding_values, Tthreshold_dict, normalization))
    log_likeli = np.mean(log_likeli)
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_df, ranges_dict)
    handles, labels = plt.gca().get_legend_handles_labels()
    SNemcee_fit_lengend = Line2D([0], [0], label=results_text, color='purple')
    handles.extend([SNemcee_fit_lengend])
    ax.legend(handles=handles, fontsize=10)
    ax.set_title('log likelihood = ' + str(round(log_likeli,3)), fontsize=18)
    ax.set_xlim(-2, 200)
    ax.set_ylim(top=43.5)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return ax


def plot_veloc_with_fit(data, sampler_df, ranges_dict, n_walkers, ax, Tthreshold_dict, normalization):
    Tthreshold = Tthreshold_dict['veloc']
    data_x = data['t_from_discovery']
    data_y = data['veloc']
    data_dy = data['dveloc']
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
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
            # y_fit_on_data_times = interp.snec_interpolator(requested[0:6], surrounding_values, models['veloc'],
            #                                                data_x_moved)
            if not isinstance(y_fit, str):
                ax.plot(x_plotting, y_fit, alpha=0.1, color='purple')
                log_likeli.append(calc_veloc_likelihood(requested, data, surrounding_values, Tthreshold_dict, normalization))
    log_likeli = np.mean(log_likeli)
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', color='k')
    ax.set_title('log likelihood = ' + str(round(log_likeli,3)), fontsize=18)
    ax.set_xlim(-2, 200)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return ax


def plot_mag_with_fit(data, sampler_df, ranges_dict, n_walkers, ax, Tthreshold_dict, normalization):
    Tthreshold = Tthreshold_dict['mag']
    filters = list(data['filter'].unique())
    data_x = data['t_from_discovery']
    y_fit = {}
    # y_fit_on_data_times = {}
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
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
                # y_fit_on_data_times[filt] = interp.snec_interpolator(requested[0:6], surrounding_values,
                #                                        models['mag'], data_x_moved, filter=filt)
                if not isinstance(y_fit, str):
                    # multiply whole graph by scaling factor
                    y_fit[filt] = y_fit[filt] * S
                    ax.plot(x_plotting, y_fit[filt], color=colors[filt], alpha=0.1)
                    log_likeli.append(calc_mag_likelihood(requested, data, surrounding_values, Tthreshold_dict, normalization))
    log_likeli = np.mean(log_likeli)
    for filt in filters:
        data_filt = data.loc[data['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_y = data_filt['abs_mag']
        data_dy = data_filt['dmag']
        ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None',
                    label=filt, color=colors[filt])
    ax.legend()
    ax.set_title('log likelihood = ' + str(round(log_likeli,3)), fontsize=18)
    ax.set_xlim(-2, 200)
    ax.invert_yaxis()
    ax.set_ylim(-14, -19)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return ax


def range_bounds(ranges_list):
    tup_list = []
    for i in range(len(ranges_list)):
        tup_list.append((np.min(ranges_list[i]), np.max(ranges_list[i])))
    return tup_list


def overlay_corner_plot(result_paths_list, output_dir, name_list, filename):
    ranges_dicts_list = []
    flat_sampler_list = []
    for result_path in result_paths_list:
        run_params_path = os.path.join(result_path, 'run_parameters.csv')
        run_params = pd.read_csv(run_params_path, index_col=0).T
        ranges_dict = import_ranges(run_params.columns.values, run_params)
        params = ranges_dict.keys()
        n_walkers = int(run_params.iloc[0]['n_walkers'])
        burn_in = int(run_params.iloc[0]['burn_in'])
        flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
        flat_sampler_noburnin = pd.read_csv(flat_sampler_path,
                                            names=params,
                                            skiprows=(burn_in - 1) * (n_walkers))
        ranges_dicts_list.append(ranges_dict)
        flat_sampler_list.append(flat_sampler_noburnin)
    labels = list(ranges_dicts_list[0].keys())
    corner_range = [1.] * len(labels)
    color_corner.overlaid_corner(
        flat_sampler_list,
        name_list,
        corner_range,
        labels,
        output_dir, filename)
    # f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
    # f_corner.savefig(os.path.join(output_dir, 'corner_plot.png'))


def get_args_from_file(result_path, ax, data_type):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
    params = ranges_dict.keys()
    initialize_empty_models_dict(ranges_dict)
    SN_name = run_params.iloc[0]['SN_name']
    n_walkers = int(run_params.iloc[0]['n_walkers'])
    burn_in = int(run_params.iloc[0]['burn_in'])
    normalization = run_params.iloc[0]['normalization']
    Tthreshold_lum = run_params.iloc[0]['Tthreshold_lum']
    Tthreshold_mag = run_params.iloc[0]['Tthreshold_mag']
    Tthreshold_veloc = run_params.iloc[0]['Tthreshold_veloc']
    Tthreshold_dict = {'lum': Tthreshold_lum, 'veloc': Tthreshold_veloc, 'mag': Tthreshold_mag}
    if data_type == 'lum':
        data = mcmc_snec.import_lum(SN_name)
    elif data_type == 'veloc':
        data = mcmc_snec.import_veloc(SN_name)
    elif data_type == 'mag':
        data = mcmc_snec.import_mag(SN_name)
    flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
    flat_sampler_noburnin = pd.read_csv(flat_sampler_path,
                                        names=params,
                                        skiprows=(burn_in - 1) * (n_walkers))
    args = (data, flat_sampler_noburnin, ranges_dict, n_walkers, ax, Tthreshold_dict, normalization)
    return args


def plot_result_fit(result_path, plot_types, ax):
    if 'lum' in plot_types:
        args = get_args_from_file(result_path, ax, 'lum')
        plot_lum_with_fit(*args)
    if 'veloc' in plot_types:
        args = get_args_from_file(result_path, ax, 'veloc')
        plot_veloc_with_fit(*args)
    if 'mag' in plot_types:
        args = get_args_from_file(result_path, ax, 'mag')
        plot_mag_with_fit(*args)
    return ax

# TODO then also need to fix the wrapper code to call stuff correctly, and also read from file names