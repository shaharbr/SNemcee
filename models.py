import numpy as np
import pandas as pd
import os


def initialize_empty_models_dict(models, ranges_dict):
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
    return models


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


def load_surrounding_models_local(models, requested, ranges_dict, fitting_type, LumTthreshold):
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
                            if LumTthreshold:
                                if models['temp'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['temp'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'temp')
    return models



