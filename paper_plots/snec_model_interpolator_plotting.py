import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rc('image', cmap='hsv')
import pandas as pd
import os

def get_edges(param_dict, param, rounding):
    return str(round(param_dict[param]['below'],rounding))+'/'\
           +str(round(param_dict[param]['above'],rounding))+' '

def snec_interpolator(requested, surrounding_values, models_dict, data_days, extend_tail=False):
    if extend_tail is not False:
        model_days = np.linspace(0, 200+extend_tail, int(1+10*(200+extend_tail)))
    else:
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
    plt.figure(figsize=(20,17))
    label_flag=True
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
                                return 'failed SN'
                            else:
                                snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir]['snec'] = \
                                    np.interp(data_days, model_days, snec_model)
                                name = 'M:' + get_edges(param_dict, 'Mzams', 1)+ \
                                       'Ni=' + get_edges(param_dict, 'Mzams', 2)+ \
                                       'E=' + get_edges(param_dict, 'Mzams', 1)+ \
                                       'R=' + get_edges(param_dict, 'Mzams', 1)+ \
                                       'K=' + get_edges(param_dict, 'Mzams', 1)+ \
                                       'Mix=' + get_edges(param_dict, 'Mzams', 1)
                                if label_flag:
                                    plt.plot(data_days,  np.log10(snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir]['snec']),
                                             label=name, color='tab:purple', alpha=0.2)
                                    label_flag = False
                                else:
                                    plt.plot(data_days,
                                             np.log10(snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir]['snec']),
                                             color='tab:purple', alpha=0.2)
                        label_flag = True
                        Mix_below = snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['below']['snec']
                        Mix_above = snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['above']['snec']
                        Mix_requested = Mix_below * param_dict['Mix']['weight_below'] + Mix_above * param_dict['Mix'][
                            'weight_above']
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['requested']['snec'] = Mix_requested

                        name = 'M:' + get_edges(param_dict, 'Mzams', 1) + \
                               'Ni=' + get_edges(param_dict, 'Mzams', 2) + \
                               'E=' + get_edges(param_dict, 'Mzams', 1) + \
                               'R=' + get_edges(param_dict, 'Mzams', 1) + \
                               'K=' + get_edges(param_dict, 'Mzams', 1) + \
                               'Mix=' + str(round(param_dict['Mix']['requested'],1))
                        if label_flag:
                            plt.plot(data_days, np.log10(Mix_requested), label=name, color='tab:blue', alpha=0.3)
                            label_flag=False
                        else:
                            plt.plot(data_days, np.log10(Mix_requested), color='tab:blue', alpha=0.3)
                    label_flag = True
                    K_below = snec_dict[Mdir][Nidir][Edir][Rdir]['below']['requested']['snec']
                    K_above = snec_dict[Mdir][Nidir][Edir][Rdir]['above']['requested']['snec']
                    K_requested = K_below * param_dict['K']['weight_below'] + K_above * param_dict['K'][
                        'weight_above']
                    snec_dict[Mdir][Nidir][Edir][Rdir]['requested']['requested']['snec'] = K_requested
                    name = 'M:' + get_edges(param_dict, 'Mzams', 1) + \
                           'Ni=' + get_edges(param_dict, 'Mzams', 2) + \
                           'E=' + get_edges(param_dict, 'Mzams', 1) + \
                           'R=' + get_edges(param_dict, 'Mzams', 1) + \
                           'K=' + str(round(param_dict['K']['requested'],1)) + ' ' + \
                           'Mix=' + str(round(param_dict['Mix']['requested'], 1))

                    if label_flag:
                        plt.plot(data_days,  np.log10(K_requested), label=name, color='tab:green', alpha=0.5)
                        label_flag=False
                    else:
                        plt.plot(data_days, np.log10(K_requested), color='tab:green', alpha=0.5)
                label_flag = True
                R_below = snec_dict[Mdir][Nidir][Edir]['below']['requested']['requested']['snec']
                R_above = snec_dict[Mdir][Nidir][Edir]['above']['requested']['requested']['snec']
                R_requested = R_below * param_dict['R']['weight_below'] + R_above * param_dict['R'][
                    'weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['requested']['requested']['snec'] = R_requested
                name = 'M:' + get_edges(param_dict, 'Mzams', 1) + \
                       'Ni=' + get_edges(param_dict, 'Mzams', 2) + \
                       'E=' + get_edges(param_dict, 'Mzams', 1) + \
                       'R=' + str(round(param_dict['R']['requested'],1)) + ' ' + \
                       'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
                       'Mix=' + str(round(param_dict['Mix']['requested'], 1))
                if label_flag:
                    plt.plot(data_days,  np.log10(R_requested), label=name, color='tab:olive', alpha=0.7)
                    label_flag=False
                else:
                    plt.plot(data_days, np.log10(R_requested), color='tab:olive', alpha=0.7)
            label_flag = True
            E_below = snec_dict[Mdir][Nidir]['below']['requested']['requested']['requested']['snec']
            E_above = snec_dict[Mdir][Nidir]['above']['requested']['requested']['requested']['snec']
            E_requested = E_below * param_dict['E']['weight_below'] + E_above * param_dict['E']['weight_above']
            snec_dict[Mdir][Nidir]['requested']['requested']['requested']['requested']['snec'] = E_requested

            name = 'M:' + get_edges(param_dict, 'Mzams', 1) + \
                   'Ni=' + get_edges(param_dict, 'Mzams', 2) + \
                   'E=' + str(round(param_dict['E']['requested'],1)) + ' ' + \
                   'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
                   'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
                   'Mix=' + str(round(param_dict['Mix']['requested'], 1))

            if label_flag:
                plt.plot(data_days,  np.log10(E_requested), label=name, color='tab:orange', alpha=0.7)
                label_flag=False
            else:
                plt.plot(data_days, np.log10(E_requested), color='tab:orange', alpha=0.7)
        label_flag = True
        Ni_below = snec_dict[Mdir]['below']['requested']['requested']['requested']['requested']['snec']
        Ni_above = snec_dict[Mdir]['above']['requested']['requested']['requested']['requested']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['requested']['requested']['requested']['requested']['snec'] = Ni_requested

        # name = 'M=' + str(round(M, 1)) + ' ' + \
        #        'Ni=' + str(round(param_dict['Ni']['below'], 1)) + ' ' + \
        #        'E=' + str(round(param_dict['E']['requested'], 1)) + ' ' + \
        #        'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
        #        'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
        #        'Mix=' + str(round(param_dict['Mix']['requested'], 1))
        # plt.plot(data_days, np.log10(Ni_below), label=name)
        #
        # name = 'M=' + str(round(M, 1)) + ' ' + \
        #        'Ni=' + str(round(param_dict['Ni']['above'], 1)) + ' ' + \
        #        'E=' + str(round(param_dict['E']['requested'], 1)) + ' ' + \
        #        'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
        #        'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
        #        'Mix=' + str(round(param_dict['Mix']['requested'], 1))
        # plt.plot(data_days, np.log10(Ni_above), label=name)

        name = 'M:' + get_edges(param_dict, 'Mzams', 1) + \
               'Ni=' + str(round(param_dict['Ni']['requested'],2)) + ' ' + \
               'E=' + str(round(param_dict['E']['requested'], 1)) + ' ' + \
               'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
               'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
               'Mix=' + str(round(param_dict['Mix']['requested'], 1))

        if label_flag:
            plt.plot(data_days,  np.log10(Ni_requested), label=name, color='tab:red', alpha=0.7)
            label_flag=False
        else:
            plt.plot(data_days, np.log10(Ni_requested), color='tab:red', alpha=0.7)
    M_below = snec_dict['below']['requested']['requested']['requested']['requested']['requested']['snec']
    M_above = snec_dict['above']['requested']['requested']['requested']['requested']['requested']['snec']
    M_requested = M_below * param_dict['Mzams']['weight_below'] + M_above * param_dict['Mzams']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec'] = M_requested


    # name = 'M=' + str(round(param_dict['Mzams']['below'], 1)) + ' ' + \
    #        'Ni=' + str(round(param_dict['Ni']['requested'], 1)) + ' ' + \
    #        'E=' + str(round(param_dict['E']['requested'], 1)) + ' ' + \
    #        'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
    #        'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
    #        'Mix=' + str(round(param_dict['Mix']['requested'], 1))
    # plt.plot(data_days, np.log10(M_below), label=name)
    #
    # name = 'M=' + str(round(param_dict['Mzams']['above'], 1)) + ' ' + \
    #        'Ni=' + str(round(param_dict['Ni']['requested'], 1)) + ' ' + \
    #        'E=' + str(round(param_dict['E']['requested'], 1)) + ' ' + \
    #        'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
    #        'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
    #        'Mix=' + str(round(param_dict['Mix']['requested'], 1))
    # plt.plot(data_days, np.log10(M_above), label=name)

    name = 'M=' + str(round(param_dict['Mzams']['requested'],1)) + '_' + \
           'Ni=' + str(round(param_dict['Ni']['requested'], 2)) + ' ' + \
           'E=' + str(round(param_dict['E']['requested'], 1)) + ' ' + \
           'R=' + str(round(param_dict['R']['requested'], 1)) + ' ' + \
           'K=' + str(round(param_dict['K']['requested'], 1)) + ' ' + \
           'Mix=' + str(round(param_dict['Mix']['requested'], 1))

    plt.plot(data_days, np.log10(M_requested), label='final interpolated model:\n'+ name, linewidth=4, color='k')


    # modelpath = os.path.join('..', 'all_lum_data', name, 'lum_observed.dat')
    # snec_model = pd.read_csv(modelpath,
    #                          names=['t_from_discovery', 'Lum'], sep=r'\s+')
    # time_col = snec_model['t_from_discovery'] / 86400  # sec to days
    # interp_days = np.linspace(0, 200, 2001)
    # snec_model = np.interp(interp_days, time_col, snec_model['Lum'])
    # plt.plot(data_days, np.log10(M_requested), label=name, linewidth=4, color='k')

    plt.title(name)
    plt.legend(fontsize=10)
    plt.ylim(41, 42.8)
    plt.xlim(20, 250)
    plt.tight_layout()
    plt.savefig('example_interpolator_plotting_'+name+'.png')
    plt.savefig('example_interpolator_plotting_' + name + '.svg')
    return snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec']


