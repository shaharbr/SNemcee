import numpy as np

data_filters = ['g', 'r', 'i', 'V', 'R', 'I']

def snec_interpolator(requested, surrounding_values, models_dict, data_days, filter=False):
    if len(data_days) < 2:
        print('no days to model')
        return 'no days to model'
    model_days = np.linspace(0, 200, 2001)
    params = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix']
    param_dict = {}
    for i in range(len(params)):
        param = params[i]
        below = surrounding_values[param][0]
        above = surrounding_values[param][1]
        if above == below:
            weight_below = 0
        else:
            weight_below = (above - requested[i]) / (above - below)
        param_dict[param] = {'requested': requested[i],
                             'below': below, 'above': above,
                             'weight_below': weight_below,
                             'weight_above': 1 - weight_below}
    # hierarchy of nested dict: 'Mzams', 'Ni', 'E', 'R', 'K', 'Mix. take note order in filename is different
    snec_dict = {'below': {}, 'above': {}, 'requested': {}}
    for Mdir in ['below', 'above', 'requested']:
        snec_dict[Mdir] = {'below': {}, 'above': {}, 'requested': {}}
        for Nidir in ['below', 'above', 'requested']:
            snec_dict[Mdir][Nidir] = {'below': {}, 'above': {}, 'requested': {}}
            for Edir in ['below', 'above', 'requested']:
                snec_dict[Mdir][Nidir][Edir] = {'below': {}, 'above': {}, 'requested': {}}
                for Rdir in ['below', 'above', 'requested']:
                    snec_dict[Mdir][Nidir][Edir][Rdir] = {'below': {}, 'above': {}, 'requested': {}}
                    for Kdir in ['below', 'above', 'requested']:
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir] = {'below': {}, 'above': {}, 'requested': {}}
                        for Mixdir in ['below', 'above', 'requested']:
                            snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir] = {}

    for Mdir in ['below', 'above']:
        for Nidir in ['below', 'above']:
            for Edir in ['below', 'above']:
                for Rdir in ['below', 'above']:
                    for Kdir in ['below', 'above']:
                        for Mixdir in ['below', 'above']:
                            M = param_dict['Mzams'][Mdir]
                            Ni = param_dict['Ni'][Nidir]
                            E = param_dict['E'][Edir]
                            R = param_dict['R'][Rdir]
                            K = param_dict['K'][Kdir]
                            Mix = param_dict['Mix'][Mixdir]
                            snec_model = models_dict[M][Ni][E][R][K][Mix]
                            if isinstance(snec_model, str):
                                print('failed SN')
                                return 'failed SN'
                            else:
                                if filter: # if it's mag
                                    snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir]['snec'] = \
                                        np.interp(data_days, model_days, snec_model[filter])
                                else:
                                    snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir]['snec'] = \
                                        np.interp(data_days, model_days, snec_model)
                        Mix_below = snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['below']['snec']
                        Mix_above = snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['above']['snec']
                        Mix_requested = Mix_below * param_dict['Mix']['weight_below'] + Mix_above * param_dict['Mix'][
                            'weight_above']
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['requested']['snec'] = Mix_requested
                    K_below = snec_dict[Mdir][Nidir][Edir][Rdir]['below']['requested']['snec']
                    K_above = snec_dict[Mdir][Nidir][Edir][Rdir]['above']['requested']['snec']
                    K_requested = K_below * param_dict['K']['weight_below'] + K_above * param_dict['K'][
                        'weight_above']
                    snec_dict[Mdir][Nidir][Edir][Rdir]['requested']['requested']['snec'] = K_requested
                R_below = snec_dict[Mdir][Nidir][Edir]['below']['requested']['requested']['snec']
                R_above = snec_dict[Mdir][Nidir][Edir]['above']['requested']['requested']['snec']
                R_requested = R_below * param_dict['R']['weight_below'] + R_above * param_dict['R'][
                    'weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['requested']['requested']['snec'] = R_requested
            E_below = snec_dict[Mdir][Nidir]['below']['requested']['requested']['requested']['snec']
            E_above = snec_dict[Mdir][Nidir]['above']['requested']['requested']['requested']['snec']
            E_requested = E_below * param_dict['E']['weight_below'] + E_above * param_dict['E']['weight_above']
            snec_dict[Mdir][Nidir]['requested']['requested']['requested']['requested']['snec'] = E_requested
        Ni_below = snec_dict[Mdir]['below']['requested']['requested']['requested']['requested']['snec']
        Ni_above = snec_dict[Mdir]['above']['requested']['requested']['requested']['requested']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['requested']['requested']['requested']['requested']['snec'] = Ni_requested
    M_below = snec_dict['below']['requested']['requested']['requested']['requested']['requested']['snec']
    M_above = snec_dict['above']['requested']['requested']['requested']['requested']['requested']['snec']
    M_requested = M_below * param_dict['Mzams']['weight_below'] + M_above * param_dict['Mzams']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec'] = M_requested
    return snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec']
