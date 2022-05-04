import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


output_dict = {'M':[], 'E':[], 'breakout':[], 'nsteps':[], 'maxtime': []}

for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]:
    for E in [0.5, 0.9, 1.3, 1.7]:
        # for Ni in [0.02]:
        #     for Mix in [2.0]:
        #         for K in [0]:
        #             for R in [0]:
        name = 'M' + str(Mzams) + \
           '_Ni' + str(0.02) + \
           '_E' + str(E) + \
           '_Mix' + str(2.0) + \
           '_R' + str(0) + \
           '_K' + str(0)
        output_dict['M'].append(Mzams)
        output_dict['E'].append(E)
        info_path = os.path.join('all_data', name, 'info.dat')
        for l in open(info_path, 'r').readlines():
            s = l.split()
            if (s[0] == 'Time') and (s[1] == 'of') and (s[2] == 'breakout'):
                time_shock_breakout = float(s[4])
                output_dict['breakout'].append(time_shock_breakout)
        log_path = os.path.join('all_data', name, name+'.txt')
        model_log = pd.read_csv(log_path, skiprows=20, usecols=[0,1], names=['nt', 'time'], sep=r'\s+')
        model_log['nt'] = pd.to_numeric(model_log['nt'], errors='coerce')
        model_log['time'] = pd.to_numeric(model_log['time'], errors='coerce')
        output_dict['nsteps'].append(model_log['nt'].iloc[-2])
        output_dict['maxtime'].append(model_log['time'].iloc[-2] / 86400)

output_df = pd.DataFrame(output_dict)
output_df.to_csv('output_model_sizes.csv')
breakout_pivot = output_df.pivot(index='M', columns='E', values='breakout')
nsteps_pivot = output_df.pivot(index='M', columns='E', values='nsteps')
maxtime_pivot = output_df.pivot(index='M', columns='E', values='maxtime')

breakout_pivot.to_csv('output_model_sizes_breakout_pivot.csv')
nsteps_pivot.to_csv('output_model_sizes_nsteps_pivot.csv')
maxtime_pivot.to_csv('output_model_sizes_maxtime_pivot.csv')



# print(output_df)
