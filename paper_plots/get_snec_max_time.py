import os
import datetime
import pandas as pd

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
home_dir = os.path.join('/home','sbracha','SNEC')
df = {'Mzams':[], 'Ni_mass':[], 'E_final':[], 'Ni_boundary':[], 'R_CSM':[], 'K_CSM':[], 'last_time':[], 'last_nt':[]}

def get_max_time(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM):
    if K_CSM == 0 or R_CSM == 0:
        K_CSM = 0
        R_CSM = 0
    name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
           + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
    outtext_path = os.path.join(home_dir, 'all_dat_combined_after_breakout', name, name+'.txt')
    outtext = pd.read_csv(outtext_path, skiprows=18, sep=r'\s+')
    last_time = outtext['time'][-1]
    last_nt = outtext['nt'][-1]
    print(outtext)
    print(last_time, last_nt)
    return last_time, last_nt



for E_final in [0.1, 0.3, 0.5]:
    for Mzams in [9.0, 10.0, 11.0]:
        for Ni_mass in [0.001, 0.02, 0.07, 0.12]:
            for Ni_boundary in [2.0, 8.0]:
                for K_CSM in [0, 10, 30, 60]:
                    for R_CSM in [0, 500, 1000, 2000]:
                           last_time, last_nt = get_max_time(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)
                           df['Mzams'].append(Mzams)
                           df['Ni_mass'].append(Ni_mass)
                           df['E_final'].append(E_final)
                           df['Ni_boundary'].append(Ni_boundary)
                           df['R_CSM'].append(R_CSM)
                           df['K_CSM'].append(K_CSM)
                           df['last_time'].append(last_time)
                           df['last_nt'].append(last_nt)


for E_final in [1.7, 1.3, 0.9, 0.5]:
    for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]:
        for Ni_mass in [0.001, 0.02, 0.07, 0.12]:
            for Ni_boundary in [2.0, 8.0]:
                for K_CSM in [0, 10, 30, 60]:
                    for R_CSM in [0, 500, 1000, 2000]:
                           last_time, last_nt = get_max_time(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)
                           df['Mzams'].append(Mzams)
                           df['Ni_mass'].append(Ni_mass)
                           df['E_final'].append(E_final)
                           df['Ni_boundary'].append(Ni_boundary)
                           df['R_CSM'].append(R_CSM)
                           df['K_CSM'].append(K_CSM)
                           df['last_time'].append(last_time)
                           df['last_nt'].append(last_nt)

df = pd.DataFrame[df]
df.DataFrame.drop_duplicates(inplace=True)
df.to_csv(os.path.join(home_dir, 'last_nt_time.csv'), index=False)


