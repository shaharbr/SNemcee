import os
import shutil, errno
from pathlib import Path


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def copy_data_to_dst(data_dir, file_name):
    file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/'+file_name) for name in list_names]
    dst_dir = [str('/home/sbracha/SNEC/'+data_dir+ '/'+ name) for name in list_names]
    dst_paths = [str('/home/sbracha/SNEC/'+data_dir+ '/'+ name + '/'+file_name) for name in list_names]
    for i in range(len(list_names)):
        dst_dir_path = Path(dst_dir[i])
        dst_file_path = Path(dst_paths[i])
        if not dst_file_path.exists():
            if not dst_dir_path.exists():
                print('making dir ' + dst_dir[i])
                os.mkdir(dst_dir[i])
            print('copying ' + file_paths[i])
            shutil.copy(file_paths[i], dst_paths[i])
        else:
            print('skipped ' + file_paths[i] + ', already exists')

list_names = [str('M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)\
                   + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM))
                           for E_final in [1.7, 1.3, 0.9, 0.5]
                           for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
                           for Ni_mass in [0.001]
                           for Ni_boundary in [2.0, 8.0]
                           for K_CSM in [0]
                           for R_CSM in [0]
                           ]

copy_data_to_dst('all_veloc_data', 'vel_Fe.dat')
copy_data_to_dst('all_lum_data', 'lum_observed.dat')
copy_data_to_dst('all_mag_data', 'magnitudes.dat')
copy_data_to_dst('all_temp_rad_data', 'T_eff.dat')
copy_data_to_dst('all_temp_rad_data', 'rad_photo.dat')
