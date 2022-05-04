import asyncio, os
import shutil, errno
import time
import re
import datetime
from add_wind_local_addedR import add_wind
import boto3
import tarfile
import compute_velocity as cv

s3 = boto3.client('s3')

project_dir = '/home/ubuntu/'
time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def progress_txt(str):
    with open(os.path.join(project_dir, time_now+'_asyncio_w_wind_output.txt'), "a") as text_file:
        text_file.write(str+'\n')

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def move_file(src_dir, dst_dir, pattern):
    for f in os.listdir(src_dir):
        if re.search(pattern, f):
            os.rename(os.path.join(src_dir, f),
                      os.path.join(dst_dir, f))

def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


async def snec(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore):
    # only enter if semaphore can be acquired
    async with semaphore:
        if K_CSM == 0 or R_CSM == 0:
            K_CSM = 0
            R_CSM = 0
        name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
               + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
        dir_name = project_dir + name + '/'
        if os.path.exists(dir_name) or os.path.exists(dir_name):
            progress_txt('skipped ' + name + ', already exists')
        elif E_final < 0.5 and Mzams > 12.0:
            progress_txt('skipped impossible SN ' + name)
        elif E_final > 0.5 and Ni_mass < 0.01:
            progress_txt('skipped unnecessary SN ' + name)
        else:
            progress_txt('start ' + name)
            snec_src = project_dir + 'SNEC/'
            copyanything(snec_src, dir_name)
            if K_CSM == 0 or R_CSM == 0:
                profile_src = project_dir + 'sukhbold_profiles/s' + str(int(Mzams)) + '/profiles/'
                copyanything(profile_src, dir_name + '/profiles/')
                time.sleep(2)
                os.chdir(os.path.abspath(dir_name))
                profile_name = 'profiles/s'+str(Mzams)+'.'
            else:
                profile_src = project_dir + 'sukhbold_profiles/s' + str(int(Mzams)) + '/profiles/'
                copyanything(profile_src, dir_name + '/profiles/')
                time.sleep(2)
                os.chdir(os.path.abspath(dir_name))
                add_wind(K_CSM, R_CSM, Mzams, dir_name)
                profile_name = 'profiles/' + 's{}_K{}_R{}.'.format(str(Mzams), str(K_CSM), str(R_CSM))

            # Read in the file
            with open('parameters', 'r') as file:
                filedata = file.read()
            # Replace the target string
            filedata = filedata.replace('profiles/xxx.', profile_name)
            filedata = filedata.replace('Ni_mass = xxx', 'Ni_mass = '+str(Ni_mass))
            filedata = filedata.replace('final_energy        = xxxd51', 'final_energy        = ' + str(E_final)+'d51')
            filedata = filedata.replace('Ni_boundary_mass = xxxd0', 'Ni_boundary_mass = ' + str(Ni_boundary) + 'd0')
            # Write the file out again
            with open('parameters', 'w') as file:
                file.write(filedata)

            cmd = './snec > '+name+'.txt'

            proc = await asyncio.create_subprocess_shell(cmd, stdout=None, stderr=None)

            await proc.communicate()
            data_src = dir_name+'Data/'
            shutil.copyfile(dir_name+name+'.txt', data_src+name+'.txt')

            cv.compute_vel(Mzams)

            dat_dir = name+'_dat'
            xg_dir = name+'_xg'
            os.mkdir(dat_dir)
            os.mkdir(xg_dir)
            move_file(data_src, xg_dir, '.*.xg')
            os.rename(data_src, dat_dir)

            make_tarfile(xg_dir + '.tar.gz', xg_dir)
            make_tarfile(dat_dir+'.tar.gz', dat_dir)
            s3.upload_file(xg_dir+'.tar.gz', 'snecbucket', xg_dir+'.tar.gz')
            s3.upload_file(dat_dir + '.tar.gz', 'snecbucket', dat_dir + '.tar.gz')
            shutil.rmtree(dir_name)
            progress_txt('end ' + name)



async def main():
    semaphore = asyncio.BoundedSemaphore(1)
    await asyncio.gather(*[snec(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore)
                           for Mzams in [12.0]
                           for E_final in [0.5]
                           for Ni_mass in [0.12]
                           for Ni_boundary in [8.0]
                           for K_CSM in [30]
                           for R_CSM in [1000]
                           ])


asyncio.run(main())
