import pandas as pd
from matplotlib import pyplot as plt
import os
import matplotlib.cm as cm


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
    resfile_CSMwith['name'] = SN_name
    resdf_CSMwith = resdf_CSMwith.append(resfile_CSMwith)
    resfile_CSMwithout = pd.read_csv(os.path.join('..', 'mcmc_results', '22_8_24', '700step',
                                                  SN_name + '_lum_csm-without_normalizedFalse_TreshLumFalse',
                                                  'final_results.csv'))
    resfile_CSMwithout['name'] = SN_name
    resdf_CSMwithout = resdf_CSMwithout.append(resfile_CSMwithout)
    paramfile = pd.read_csv(os.path.join('..', 'mcmc_results', '22_8_24', '700step', SN_name+'_lum_csm-with_normalizedFalse_TreshLumFalse', 'run_parameters.csv'), index_col=0).T


# resdf = resdf.sort_values(by='name')
# resdf_noCSM = resdf_noCSM.sort_values(by='name')
resdf_CSMwithout = resdf_CSMwithout.add_suffix('_noCSM')
resdf_CSMwithout = resdf_CSMwithout.rename(columns={'name_noCSM':'name'})



resdf = resdf_CSMwith.merge(resdf_CSMwithout, on='name')

def plot_csm_scatter(param, param_range, ax):
    y_list = list(resdf[param+'_noCSM'])
    x_list = list(resdf[param])
    y_lower = list(resdf[param+'_lower_noCSM'])
    x_lower = list(resdf[param+'_lower'])
    y_upper = list(resdf[param+'_upper_noCSM'])
    x_upper = list(resdf[param+'_upper'])
    names = [name[4:] for name in list(resdf['name'])]
    K = list(resdf['K'])
    R = list(resdf['R'])
    ax.plot(param_range, param_range, color='gray')
    ax.errorbar(x_list, y_list,
                xerr=[x_lower, x_upper], yerr=[y_lower, y_upper],
                linestyle='None', color='black')
    axes = []
    K_color_norm = 50
    K_examples = [10, 30, 60]
    R_marker_norm = 80
    R_examples = [100, 500, 1000]
    for i, txt in enumerate(names):
        x1 = min(x_list)
        x2 = max(x_list)
        xspacing = (x2 - x1) / 50
        l, = ax.plot(x_list[i], y_list[i], label=names[i], marker='o',
                     markersize=(R[i]+50)/R_marker_norm, color=cm.viridis((70-K[i])/K_color_norm),linestyle = 'None')
        ax.annotate(txt, (x_list[i]+xspacing, y_list[i]+xspacing))
        axes.append(l)
    ax.set_ylabel(param+' without CSM')
    ax.set_xlabel(param+' with CSM')
    ax.legend()
    gll, = ax.plot([],[], markersize=(R_examples[0]+50)/R_marker_norm, marker='o', color=cm.viridis((70-K_examples[0])/K_color_norm),linestyle = 'None')
    gl, = ax.plot([],[], markersize=(R_examples[1]+50)/R_marker_norm, marker='o', color=cm.viridis((70-K_examples[1])/K_color_norm),linestyle = 'None')
    ga, = ax.plot([],[], markersize=(R_examples[2]+50)/R_marker_norm, marker='o', color=cm.viridis((70-K_examples[2])/K_color_norm),linestyle = 'None')
    # legend1 = ax.legend(axes, names, loc='upper left', fontsize=12)
    # ax.add_artist(legend1)
    ax.legend((gll,gl,ga),
           ('K='+str(K_examples[0])+' R='+str(R_examples[0]),
            'K='+str(K_examples[1])+' R='+str(R_examples[1]),
            'K='+str(K_examples[2])+' R='+str(R_examples[2])),
           scatterpoints=1,
           loc='lower right',
           ncol=1,
           fontsize=12)
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


fig, axs = plt.subplots(3, 2, figsize=(12,18))

params = ['E', 'Mzams', 'Ni', 'Mix', 'S', 'T']

rows = 3
columns = 2
for r in range(rows):
    for c in range(columns):
        i = (r+1) * (c+1) - 1
        plot_csm_scatter(params[i], param_ranges[params[i]], axs[r,c])

plt.tight_layout()
plt.savefig(os.path.join('figures', 'csm_effect_scatter.png'))
plt.savefig(os.path.join('figures', 'csm_effect_scatter.pdf'))

#TODO compare one step vs two step with/wo prior carryover for a given SN, and show traces

plt.show()