import os
os.environ['PYSYN_CDBS'] = os.path.join('..', 'cdbs')
import pandas as pd
import numpy as np
import pysynphot
from datetime import datetime


home_dir = os.path.join('/home','sbracha','SNEC')
time_now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
solar_rad = 6.957 * (10 ** 10)  # cm

times = list(np.arange(0.0, 2.0, 0.1) * 86400) \
        + list(np.arange(2.0, 10.0, 0.5) * 86400) \
        + list(np.arange(10.0, 100.0, 5.0) * 86400) \
        + list(np.arange(100.0, 150.0, 0.5) * 86400) \
        + list(np.arange(150.0, 250.0, 5.0) * 86400)

num_timepoints = len(times)

names = [str('M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)
                 + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM))
                           for K_CSM in [0]
                           for E_final in [0.5]
                           for Ni_mass in [0.001]
                           for Mzams in [9.0]
                           for R_CSM in [0]
                           for Ni_boundary in [2.0]
             ]

filters = {'u': 'sdss,u', 'g': 'sdss,g', 'r': 'sdss,r', 'i': 'sdss,i', 'z': 'sdss,z',
           'U': 'landolt,u', 'B': 'landolt,b', 'V': 'landolt,v', 'R': 'landolt,r', 'I': 'landolt,i'}


def progress_txt(str):
    with open(os.path.join(home_dir, time_now+'_lum_to_mag.txt'), "a") as text_file:
        text_file.write(str+'\n')


def lum_to_mag(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM):
    if K_CSM == 0 or R_CSM == 0:
        K_CSM = 0
        R_CSM = 0
    name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
           + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)

    dst_dir = os.path.join(home_dir, 'all_dat_combined_after_breakout', name)
    dst_filepath = os.path.join(dst_dir, 'magnitudes_pys.dat')
    if os.path.exists(dst_filepath):
        progress_txt('skipped ' + name + ', pys mag file already exists')
    else:
        mag_dict = {'t_from_discovery': [], 'Teff': [], 'PTF_R_AB': [], 'R^2': []}
        for filtername in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
            mag_dict[filtername] = []

        temp_table = pd.read_csv(os.path.join('all_dat_combined_after_breakout', name, 'T_eff.dat'))
        rad_table = pd.read_csv(os.path.join('all_dat_combined_after_breakout', name, 'rad_photo.dat'))
        for time in times:
            temp = np.interp(time, temp_table['t_from_discovery'], temp_table['temp'])
            if temp > 500:
                mag_dict['t_from_discovery'].append(time)
                mag_dict['PTF_R_AB'].append(0)
                mag_dict['Teff'].append(temp)
                rad = np.interp(time, rad_table['t_from_discovery'], rad_table['rad']) / solar_rad # in solar radii
                mag_dict['R^2'].append(rad ** 2)
                bb = pysynphot.BlackBody(temp)  # For 1 Solar Radius at 1 kpc
                for filt in filters.keys():
                    filter_name = filters[filt]
                    bp = pysynphot.ObsBandpass(filter_name)
                    obs = pysynphot.Observation(bb, bp)
                    if ('sdss' in filter_name) or filter_name == 'o':
                        magsystem = 'ABMag'
                    else:
                        magsystem = 'VegaMag'
                    #     TODO plot the radius and temp over time for the models, to understand why similar lum lead to such diff mags
                    expected_mag = obs.effstim(magsystem) - 2.5*np.log10((rad**2)*((1000.0/10)**2)) # divide by 10 to get absolute magnitude (distance of 10 parsec)
                    mag_dict[filt].append(expected_mag)
            else:
                progress_txt(name, '- temp too low: ', str(temp), 'time: ', str(time))
                break
                # mag_dict[filt].append(obs.effstim(magsystem)) # Rescaling from the default (1 solar radius at 1000 pc)
                # TODO sort out whats going on here - why is this the formula, why is dist_mod and redshif removed, what the correct thing to do
        mag_df = pd.DataFrame(mag_dict)
        # plt.plot(mag_df['t_from_discovery']/86400, mag_df['V'] - 36.01, marker='o')
        # plt.gca().invert_yaxis()
        # plt.show()
        mag_df.to_csv(os.path.join(dst_dir, 'magnitudes_pys.dat'),
                      header=False, index=False)



for E_final in [0.1, 0.3, 0.5]:
    for Mzams in [9.0, 10.0, 11.0]:
        for Ni_mass in [0.001, 0.02, 0.07, 0.12]:
            for Ni_boundary in [2.0, 8.0]:
                for K_CSM in [0, 10, 30, 60]:
                    for R_CSM in [0, 500, 1000, 2000]:
                           lum_to_mag(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)


for E_final in [1.7, 1.3, 0.9, 0.5]:
    for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]:
        for Ni_mass in [0.001, 0.02, 0.07, 0.12]:
            for Ni_boundary in [2.0, 8.0]:
                for K_CSM in [0, 10, 30, 60]:
                    for R_CSM in [0, 500, 1000, 2000]:
                           lum_to_mag(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)