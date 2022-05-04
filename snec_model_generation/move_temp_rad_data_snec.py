import os
import shutil, errno


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


list_names = [str('M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)\
                   + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM))
                           for K_CSM in [60]
                           for E_final in [1.7, 1.3, 0.9, 0.5]
                           for Ni_mass in [0.12, 0.02, 0.07]
                           for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
                           for R_CSM in [500, 1000]
                           for Ni_boundary in [2.0, 8.0]
              ]


# SNEC_dir = os.path.join('home', 'sbracha', 'SNEC')
# dst_dir = [os.path.join(SNEC_dir, 'all_temp_rad_data', name) for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_temp_rad_data/' + name) for name in list_names]

# T_file_paths = [os.path.join(SNEC_dir, 'all_data', name, 'T_eff.dat') for name in list_names]
T_file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/T_eff.dat') for name in list_names]
# T_dst_paths = [os.path.join(SNEC_dir, 'all_temp_rad_data', name, 'T_eff.dat') for name in list_names]
T_dst_paths = [str('/home/sbracha/SNEC/all_temp_rad_data/' + name + '/T_eff.dat') for name in list_names]

# R_file_paths = [os.path.join(SNEC_dir, 'all_data', name, 'rad_photo.dat') for name in list_names]
R_file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/rad_photo.dat') for name in list_names]
# R_dst_paths = [os.path.join(SNEC_dir, 'all_temp_rad_data', name, 'rad_photo.dat') for name in list_names]
R_dst_paths = [str('/home/sbracha/SNEC/all_temp_rad_data/' + name + '/rad_photo.dat') for name in list_names]


for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(T_file_paths[i], T_dst_paths[i])
    shutil.copy(R_file_paths[i], R_dst_paths[i])
