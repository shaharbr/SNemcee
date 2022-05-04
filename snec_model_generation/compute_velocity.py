import os
import pandas as pd
from pylab import *
import xgraph as xg
import math

msol = 1.98892e33
rsol = 6.955e10
clite = 2.998e10
m_e = 9.1094e-28
e_e = 4.803e-10
amu = 1.6605e-24
Fe_mass = 55.845*amu
ff = 0.023
lambda_0 = 5169e-8


def compute_vel(Mzams, dat_dir, xg_dir, home_dir):
    rho = xg.parsefile(os.path.join(xg_dir, 'rho.xg'))
    temp = xg.parsefile(os.path.join(xg_dir, 'temp.xg'))
    vel = xg.parsefile(os.path.join(xg_dir, 'vel.xg'))

    profile_kepler = os.path.join(home_dir, 'sukhbold_profiles_kepler', 's'+str(Mzams)+'_presn')
    mass_profile, Fe_mass_fraction = np.loadtxt(profile_kepler,usecols=(1,30),unpack=True,skiprows=2)

    n_Fe = []
    for l in range(rho.nframes):
        n_Fe_x = []
        for k in range(len(rho.frame(l).data_x)):
            mass_fraction_current = np.interp(rho.frame(l).data_x[k],mass_profile,Fe_mass_fraction)
            n_Fe_x.append(rho.frame(l).data_y[k]*mass_fraction_current/Fe_mass)
        n_Fe.append(n_Fe_x)

    T_tab = []
    rho_tab = []

    eta_file = os.path.join(home_dir, 'FeII_5169_eta.dat')

    for l in open(eta_file, 'r').readlines():
        s = l.split()
        if s[0] == '1e-16':
            T_tab.append(float(s[1]))
        if s[1] == '2000':
            rho_tab.append(float(s[0]))

    eta_tab = []
    eta_tab_x = []

    k = 0

    for l in open(eta_file, 'r').readlines():
        s = l.split()

        eta_tab_x.append(float(s[2]))
        if k == len(T_tab) - 1:
            eta_tab.append(eta_tab_x)
            eta_tab_x = []
            k = 0
        else:
            k = k + 1

    found_time_break = 0
    info_path = os.path.join(dat_dir, 'info.dat')
    for l in open(info_path, 'r').readlines():
        s = l.split()
        if (s[0] == 'Time') and (s[1] == 'of') and (s[2] == 'breakout'):
            time_shock_breakout = float(s[4])
            found_time_break = 1
        if found_time_break == 0:
            time_shock_breakout = 0.0

    index_phot_path = os.path.join(dat_dir, 'index_photo.dat')

    time_ind_ph, ind_ph = np.loadtxt(index_phot_path,usecols=(0,1),unpack=True)

    vel_tau_sob = []
    time_tau_sob = []

    for l in range(vel.nframes):

        time_tau_sob.append((vel.frame(l).time-time_shock_breakout)/24.0/3600.0)

        index_photo_curr = np.interp(vel.frame(l).time,time_ind_ph,ind_ph)
        found_tau_sob = 0
        vel_tau_sob_x = 0.0

        for i in range(int(index_photo_curr),len(vel.frame(l).data_x)):
            temp_current = temp.frame(l).data_y[i]
            rho_current = rho.frame(l).data_y[i]
            eta_temp = []
            for k in range(len(eta_tab)):
                eta_temp.append(np.interp(temp_current,T_tab,eta_tab[k]))
            eta_current = np.interp(rho_current,rho_tab,eta_temp)
            tau_sob_curr = math.pi*e_e*e_e/m_e/clite*n_Fe[l][i]*eta_current*ff*(rho.frame(l).time-time_shock_breakout)*lambda_0
            if tau_sob_curr < 1.0 and found_tau_sob == 0:
                vel_tau_sob_x = vel.frame(l).data_y[i]*1e-5
                found_tau_sob = 1

        vel_tau_sob.append(vel_tau_sob_x)

    sob1 = pd.DataFrame({'t_from_discovery': np.array(time_tau_sob), 'veloc': vel_tau_sob})
    sob1.to_csv(os.path.join(dat_dir, 'vel_Fe.dat'))

