import pandas as pd
from matplotlib import pyplot as plt
import os
import seaborn as sns
import matplotlib.cm as cm
import numpy as np


resdf_CSMwith = pd.DataFrame({'name':[],
                      'Mzams':[],'Mzams_lower':[],'Mzams_upper':[],
                      'Ni':[],'Ni_lower':[],'Ni_upper':[],
                      'E':[],'E_lower':[],'E_upper':[],
                      'R':[],'R_lower':[],'R_upper':[],
                      'K':[],'K_lower':[],'K_upper':[],
                      'Mix':[],'Mix_lower':[],'Mix_upper':[],
                      'T':[],'T_lower':[],'T_upper':[],
                      'S':[],'S_lower':[],'S_upper':[]})

resdf_CSMwithout = pd.DataFrame({'name':[],
                      'Mzams':[],'Mzams_lower':[],'Mzams_upper':[],
                      'Ni':[],'Ni_lower':[],'Ni_upper':[],
                      'E':[],'E_lower':[],'E_upper':[],
                      'Mix':[],'Mix_lower':[],'Mix_upper':[],
                      'T':[],'T_lower':[],'T_upper':[],
                      'S':[],'S_lower':[],'S_upper':[]})

SN_names = ['SN2004a', 'SN2005cs', 'SN2008bk', 'SN2012aw', 'SN2012ec', 'SN2017eaw', 'SN2018aoq', 'SN2020bij']

for SN_name in SN_names:
    print(SN_name + '_lum_csm-with_normalizedFalse_TreshLumFalse')
    resfile_CSMwith = pd.read_csv(os.path.join('..', 'mcmc_results', '22_8_24', '700step',
                                                  SN_name + '_lum_csm-with_normalizedFalse_TreshLumFalse',
                                                  'final_results.csv'))
    print(resfile_CSMwith)
    resfile_CSMwith['name'] = SN_name
    resdf_CSMwith.append(resfile_CSMwith)

    resfile_CSMwithout = pd.read_csv(os.path.join('..', 'mcmc_results', '22_8_24', '700step',
                                                  SN_name + '_lum_csm-without_normalizedFalse_TreshLumFalse',
                                                  'final_results.csv'))
    resfile_CSMwithout['name'] = SN_name
    resdf_CSMwithout.append(resfile_CSMwithout)
    # print(resdf_CSMwithout)

    paramfile = pd.read_csv(os.path.join('..', 'mcmc_results', '22_8_24', '700step', SN_name+'_lum_csm-with_normalizedFalse_TreshLumFalse', 'run_parameters.csv'), index_col=0).T


# resdf = resdf.sort_values(by='name')
# resdf_noCSM = resdf_noCSM.sort_values(by='name')
resdf_CSMwithout = resdf_CSMwithout.add_suffix('_noCSM')
resdf_CSMwithout = resdf_CSMwithout.rename(columns={'name_noCSM':'name'})



resdf = resdf_CSMwith.merge(resdf_CSMwithout, on='name')

def plot_csm_scatter(param, csm_param, param_range, ax):
    y_list = list(resdf[param+'_noCSM'])
    x_list = list(resdf[param])
    y_lower = list(resdf[param+'_lower_noCSM'])
    x_lower = list(resdf[param+'_lower'])
    y_upper = list(resdf[param+'_upper_noCSM'])
    x_upper = list(resdf[param+'_upper'])
    names = list(resdf['name'])
    CSM = list(resdf[csm_param])
    ax.plot(param_range, param_range, color='gray')
    ax.errorbar(x_list, y_list,
                xerr=[x_lower, x_upper], yerr=[y_lower, y_upper],
                linestyle='None', color='black')
    axes = []
    if csm_param == 'K':
        marker_norm = 4
        color_norm = 50
        param_examples = [10, 30, 60]
    elif csm_param == 'R':
        marker_norm = 60
        color_norm = 700
        param_examples = [100, 500, 1000]
    for i in range(len(names)):
        l, = ax.plot(x_list[i], y_list[i], label=names[i], marker='o',
                     markersize=CSM[i]/marker_norm, color=cm.viridis(CSM[i]/color_norm),linestyle = 'None')
        axes.append(l)
    ax.set_ylabel(param+' without CSM')
    ax.set_xlabel(param+' with CSM')
    ax.legend()
    gll, = ax.plot([],[], markersize=param_examples[0]/marker_norm, marker='o', color=cm.viridis(param_examples[0]/color_norm),linestyle = 'None')
    gl, = ax.plot([],[], markersize=param_examples[1]/marker_norm, marker='o', color=cm.viridis(param_examples[1]/color_norm),linestyle = 'None')
    ga, = ax.plot([],[], markersize=param_examples[2]/marker_norm, marker='o', color=cm.viridis(param_examples[2]/color_norm),linestyle = 'None')
    legend1 = ax.legend(axes, names, loc='upper left', fontsize=12)
    ax.legend((gll,gl,ga),
           (str(param_examples[0]), str(param_examples[1]), str(param_examples[2])),
           scatterpoints=1,
           title=csm_param,
           loc='lower left',
           ncol=1,
           fontsize=12)
    ax.add_artist(legend1)
    # plt.savefig(os.path.join('figures', 'csm_effect_scatter_'+param+'_'+csm_param+'.png'))
    # plt.savefig(os.path.join('figures', 'csm_effect_scatter_' + param + '_' + csm_param + '.svg'))




Mzams_range = [9.0, 10.0, 11.0, 13.0, 15.0, 17.0]  # Progenitor ZAMS mass (solar mass)
Ni_range = [0.001, 0.02, 0.07, 0.12]  # Ni mass (solar mass)
E_final_range = [0.1, 0.3, 0.5, 0.9, 1.3, 1.7]  # Explosion energy (foe=10^51 erg)
Mix_range = [2.0, 8.0]  # Ni mixing (solar mass)
R_range = [0, 500, 1000]
K_range = [0, 10, 30, 60]
T_range = [-5, 2]  # time shift from assumed day of explosion (days)
param_ranges = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                    'R': R_range, 'K': K_range, 'Mix': Mix_range, 'S': [], 'T': T_range}


fig, axs = plt.subplots(6, 2, figsize=(10,25))

params = ['E', 'Mzams', 'Ni', 'Mix', 'S', 'T']
csm_params = ['K', 'R']

for i in range(6):
    for j in range(2):
        plot_csm_scatter(params[i], csm_params[j], param_ranges[params[i]], axs[i,j])

plt.tight_layout()
plt.savefig(os.path.join('figures', 'csm_effect_scatter.png'))
plt.savefig(os.path.join('figures', 'csm_effect_scatter.svg'))



#TODO check on the new SNEC runs?

#TODO compare one step vs two step with/wo prior carryover for a given SN, and show traces

plt.show()