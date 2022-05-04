import os
import pandas as pd
from datetime import datetime
import shutil, errno


home_dir = os.path.join('/home','sbracha','SNEC')
time_now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

vel_col = ['unnamed', 't_from_discovery', 'vel']
lum_col = ['t_from_discovery', 'Lum']
mag_col = ['t_from_discovery', 'Teff', 'PTF_R_AB', 'R^2', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']
temp_col = ['t_from_discovery', 'temp']
rad_col = ['t_from_discovery', 'rad']


def progress_txt(str):
    with open(os.path.join(home_dir, time_now+'_remove_time_before_breakout.txt'), "a") as text_file:
        text_file.write(str+'\n')

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def remove_before_breakout(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM):
    if K_CSM == 0 or R_CSM == 0:
        K_CSM = 0
        R_CSM = 0
    name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
           + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)

    dst_dir = os.path.join(home_dir, 'all_dat_combined_after_breakout', name)
    if os.path.exists(dst_dir):
        if R_CSM != 0:
            progress_txt('skipped ' + name + ', dir already exists')
    else:
        os.mkdir(dst_dir)
    info_path = os.path.join(home_dir, 'all_dat_combined', name, 'info.dat')
    found_time_break = 0
    for l in open(info_path, 'r').readlines():
        s = l.split()
        if (s[0] == 'Time') and (s[1] == 'of') and (s[2] == 'breakout'):
            time_shock_breakout = float(s[4])
            found_time_break = 1
        if found_time_break == 0:
            time_shock_breakout = 0.0

    for filename, column_names in zip(['vel_Fe.dat', 'lum_observed.dat', 'magnitudes.dat', 'T_eff.dat', 'rad_photo.dat'],
                                      [vel_col, lum_col, mag_col, temp_col, rad_col]):
        og_dir = os.path.join(home_dir, 'all_dat_combined', name)
        dst_dir = os.path.join(home_dir, 'all_dat_combined_after_breakout', name)


        filepath = os.path.join(og_dir, filename)
        if filename == 'vel_Fe.dat':
            df = pd.read_csv(filepath, names=column_names)
            df.drop(columns=['unnamed'], inplace=True)
            df.drop(index=[0], inplace=True)
            df['t_from_discovery'] = df['t_from_discovery'].astype(float)
            df['t_from_discovery'] = df['t_from_discovery'] * 86400  # sec to days
        else:
            df = pd.read_csv(filepath, names=column_names, sep=r'\s+')
            df['t_from_discovery'] = df['t_from_discovery'] - time_shock_breakout
        df = df.loc[df['t_from_discovery'] > 0]
        df.to_csv(os.path.join(home_dir, 'all_dat_combined_after_breakout', name, filename), index=False)

    copyanything(os.path.join(og_dir, name + '.txt'), dst_dir)




for E_final in [0.1, 0.3]:
    for Mzams in [12.0, 13.0, 15.0, 17.0]:
        for Ni_mass in [0.07]:
            for Ni_boundary in [2.0]:
                for K_CSM in [0]:
                    for R_CSM in [0]:
                           remove_before_breakout(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)
