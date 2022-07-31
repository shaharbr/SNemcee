import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import os
import mcmc_snec
import re
import colorful_corner_plot as color_corner


T_thresh = 10 ** 3.5
extend_tail = False
filters = ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']
colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}


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


def open_reshape_3d_array(output_dir, type, step):
    array_2d = np.genfromtxt(os.path.join(output_dir, type+ '_models_step' + str(step) + '.txt'))
    shape_2d = array_2d.shape
    array_3d = array_2d.reshape((shape_2d[0], 2, shape_2d[1]/2))
    return array_3d


def plot_lum_with_fit(data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumTthreshold):
    data = data_dict['lum']
    data_x = data['t_from_discovery']
    data_y = data['Lum']
    dy0 = data['dLum0']
    dy1 = data['dLum1']
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            data_x_moved = data_x - T
            if LumTthreshold:
                max_x, temp_fit = mcmc_snec.temp_thresh_cutoff(requested[0:6], ranges_dict, data_x_moved)
                data = data.loc[data_x_moved <= max_x]
                x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            else:
                x_plotting = np.linspace(-T, 200-T, 2001)
            y_fit = mcmc_snec.interp_yfit(requested, ranges_dict, 'lum', x_plotting)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                # y_fit_on_data_times = y_fit_on_data_times * S
                ax.plot(x_plotting, np.log10(y_fit), alpha=0.1, color='purple')
                log_likeli.append(mcmc_snec.log_likelihood(requested, data_dict, ranges_dict, 'lum', LumTthreshold, normalization))
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


def plot_veloc_with_fit(data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumTthreshold):
    data = data_dict['veloc']
    data_x = data['t_from_discovery']
    data_y = data['veloc']
    data_dy = data['dveloc']
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            data_x_moved = data_x - T
            # always truncate by temp thresh for veloc
            max_x, temp_fit = mcmc_snec.temp_thresh_cutoff(requested[0:6], ranges_dict, data_x_moved)
            data = data.loc[data_x_moved <= max_x]
            x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            x_plotting_moved = x_plotting - T
            y_fit = mcmc_snec.interp_yfit(requested, ranges_dict, 'veloc', x_plotting)
            if not isinstance(y_fit, str):
                ax.plot(x_plotting, y_fit, alpha=0.1, color='purple')
                log_likeli.append(
                    mcmc_snec.log_likelihood(requested, data_dict, ranges_dict, 'veloc', LumTthreshold, normalization))
    log_likeli = np.mean(log_likeli)
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', color='k')
    ax.set_title('log likelihood = ' + str(round(log_likeli,3)), fontsize=18)
    ax.set_xlim(-2, 200)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return ax


def plot_mag_with_fit(data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumTthreshold):
    data = data_dict['mag']
    filters = list(data['filter'].unique())
    data_x = data['t_from_discovery']
    y_fit = {}
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            data_x_moved = data_x - T
            # always truncate by temp thresh for mag
            max_x, temp_fit = mcmc_snec.temp_thresh_cutoff(requested[0:6], ranges_dict, data_x_moved)
            data = data.loc[data_x_moved <= max_x]
            x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            for filt in filters:
                y_fit[filt] = mcmc_snec.interp_yfit(requested, ranges_dict, 'mag', x_plotting, filter=filt)
                if not isinstance(y_fit, str):
                    # multiply whole graph by scaling factor
                    y_fit[filt] = y_fit[filt] -2.5*np.log10(S)
                    ax.plot(x_plotting, y_fit[filt], color=colors[filt], alpha=0.1)

                    log_likeli.append(
                        mcmc_snec.log_likelihood(requested, data_dict, ranges_dict, 'mag', LumTthreshold,
                                                 normalization))
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


def string_to_bool(mystr):
    if mystr == 'False':
        return False
    else:
        return True

# TODO adapt these two below to work in plot_snec_fits (it was in mcmc_snec)

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

def chain_plots(result_path, output_dir, first_stage_steps=None):
    # (result_paths_list, output_dir, name_list, filename):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
    params = ranges_dict.keys()
    n_walkers = int(run_params.iloc[0]['n_walkers'])
    n_steps = int(run_params.iloc[0]['n_steps'])
    burn_in = int(run_params.iloc[0]['burn_in'])
    flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
    sampler_array = np.genfromtxt(flat_sampler_path, delimiter=',')
    sampler_chain = sampler_array.reshape((n_walkers, n_steps, len(params)))



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


def get_args_from_file(result_path, ax, data_type):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
    params = ranges_dict.keys()
    SN_name = run_params.iloc[0]['SN_name']
    # TODO shoud this be n_steps or n_steps-1?
    # step = run_params.iloc[0]['n_steps']
    n_walkers = int(run_params.iloc[0]['n_walkers'])
    burn_in = int(run_params.iloc[0]['burn_in'])
    normalization = run_params.iloc[0]['normalization']
    LumTthreshold = string_to_bool(run_params.iloc[0]['Tthreshold_lum'])
    data_dict = mcmc_snec.load_SN_data(data_type, SN_name)
    flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
    sampler_df = pd.read_csv(flat_sampler_path,
                                        names=params,
                                        skiprows=(burn_in - 1) * (n_walkers))
    args = (data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumTthreshold)
    return args


def plot_result_fit(result_path, plot_types, ax):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
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