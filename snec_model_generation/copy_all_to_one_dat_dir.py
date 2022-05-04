import os
import shutil, errno
import datetime

project_dir = '/home/sbracha/SNEC/'
time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def progress_txt(str):
    with open(os.path.join(project_dir, time_now+'_copy_all_to_one_dat_dir.txt'), "a") as text_file:
        text_file.write(str+'\n')


home_dir = os.path.join('/home','sbracha','SNEC')

def transfer_all_to_data_dir(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM):
    if K_CSM == 0 or R_CSM == 0:
        K_CSM = 0
        R_CSM = 0
    name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
           + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
    dst_dir = os.path.join(home_dir, 'all_dat_combined', name)
    if os.path.exists(dst_dir):
        progress_txt('skipped ' + name + ', dir already exists')
    else:
        os.mkdir(dst_dir)
    for filename in ['vel_Fe.dat', 'lum_observed.dat', 'magnitudes.dat', 'T_eff.dat', 'rad_photo.dat', 'info.dat', name+'.txt']:
        og_filepath = os.path.join(home_dir, 'all_data', name, filename)
        dst_filepath = os.path.join(dst_dir, filename)
        if os.path.exists(dst_filepath):
            progress_txt('skipped ' + name + ', file already copied there')
        else:
            copyanything(og_filepath, dst_dir)




for E_final in [0.1, 0.3, 0.5]:
    for Mzams in [9.0, 10.0, 11.0]:
        for Ni_mass in [0.001, 0.02, 0.07, 0.12]:
            for Ni_boundary in [2.0, 8.0]:
                for K_CSM in [0, 10, 30, 60]:
                    for R_CSM in [0, 500, 1000, 2000]:
                           transfer_all_to_data_dir(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)


for E_final in [1.7, 1.3, 0.9, 0.5]:
    for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]:
        for Ni_mass in [0.001, 0.02, 0.07, 0.12]:
            for Ni_boundary in [2.0, 8.0]:
                for K_CSM in [0, 10, 30, 60]:
                    for R_CSM in [0, 500, 1000, 2000]:
                           transfer_all_to_data_dir(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM)

