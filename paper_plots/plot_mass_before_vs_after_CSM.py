import pandas as pd
from matplotlib import pyplot as plt
import os
import seaborn as sns
import matplotlib.cm as cm
import numpy as np



list_results = ['2021-07-14_05-55-48_lum_SN2004a',
                '2021-07-14_07-17-52_lum_SN2012ec',
                '2021-07-14_07-03-07_lum_noCSM_SN2017eaw',
                '2021-07-14_07-01-35_lum_SN2012aw',
                '2021-07-14_06-51-58_lum_noCSM_SN2012ec',
                '2021-07-14_06-41-38_lum_noCSM_SN2012aw',
                '2021-07-14_06-11-00_lum_SN2004et_d5.9',
                '2021-07-14_06-08-05_lum_noCSM_SN2004et_d5.9',
                '2021-07-14_05-58-39_lum_noCSM_SN2004a',
                '2021-07-14_07-34-21_lum_SN2017eaw']

resdf = pd.DataFrame({'name':[],
                      'Mzams':[],'Mzams_lower':[],'Mzams_upper':[],
                      'Ni':[],'Ni_lower':[],'Ni_upper':[],
                      'E':[],'E_lower':[],'E_upper':[],
                      'R':[],'R_lower':[],'R_upper':[],
                      'K':[],'K_lower':[],'K_upper':[],
                      'Mix':[],'Mix_lower':[],'Mix_upper':[],
                      'T':[],'T_lower':[],'T_upper':[],
                      'S':[],'S_lower':[],'S_upper':[]})

resdf_noCSM = pd.DataFrame({'name':[],
                      'Mzams':[],'Mzams_lower':[],'Mzams_upper':[],
                      'Ni':[],'Ni_lower':[],'Ni_upper':[],
                      'E':[],'E_lower':[],'E_upper':[],
                      'Mix':[],'Mix_lower':[],'Mix_upper':[],
                      'T':[],'T_lower':[],'T_upper':[],
                      'S':[],'S_lower':[],'S_upper':[]})

for resdir in list_results:
    resfile = pd.read_csv(os.path.join('..', '2021-07-14', resdir, 'final_results.csv'))
    paramfile = pd.read_csv(os.path.join('..', '2021-07-14', resdir, 'run_parameters.csv'), index_col=0).T
    resfile['name'] = str(paramfile.iloc[0]['SN_name'])
    if 'noCSM' in resdir:
        resdf_noCSM = resdf_noCSM.append(resfile)
    else:
        resdf = resdf.append(resfile)

# resdf = resdf.sort_values(by='name')
# resdf_noCSM = resdf_noCSM.sort_values(by='name')
resdf_noCSM = resdf_noCSM.add_suffix('_noCSM')
resdf_noCSM = resdf_noCSM.rename(columns={'name_noCSM':'name'})

print(resdf)
print(resdf_noCSM)

resdf = resdf.merge(resdf_noCSM, on='name')

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
        param_examples = [10, 30, 50]
    elif csm_param == 'R':
        marker_norm = 60
        color_norm = 700
        param_examples = [200, 500, 800]
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

param_ranges = {'E':[0.4, 1.5], 'Mzams':[9, 17],
                'Ni':[0.02, 0.12], 'Mix':[2, 8],
                'S':[0.7, 1.2], 'T':[0,1.2]}



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