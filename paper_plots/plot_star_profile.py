import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt


'''
mass: gram
radius: cm
temp: K
rho: gr / cm^3
vel: cm / s

'''

s15 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's15.0.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

s15_K60_R1000 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's15.0_K60_R1000.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

m_sol = 1.989 * (10 ** 33)  # gram
r_sol = 6.9634 * (10 ** 10)  # cm

# plt.figure(figsize=(3,2))
# plt.plot(s15['radius'] / r_sol, s15['rho'], color='k')
# plt.ylabel('log rho')
# plt.xlabel('radius (in R_sol)')
# plt.title('s15_kepler')
# plt.yscale('log')
#
# plt.figure(figsize=(3,2))
# plt.plot(s15['mass'] / m_sol, s15['rho'], color='k')
# plt.ylabel('log rho')
# plt.xlabel('mass (in M_sol)')
# plt.title('s15_kepler')
# plt.yscale('log')


# print
K = 60
R = 1000
r_vec = s15_K60_R1000['radius']

csm_rho = np.ones(len(r_vec)) * K *(10**17) / (r_vec**2)
# csm_rho_integral = 4 * np.pi * K*(10**17) * R # volumetric spheric integral
# print(csm_rho_integral / m_sol)


# plt.figure(figsize=(3,2))
# plt.plot(s15_K60_R1000['radius'] / r_sol, csm_rho)
# plt.ylabel('log rho')
# plt.xlabel('radius (in R_sol)')
# plt.title('CSM')
# plt.yscale('log')
#
# plt.figure(figsize=(3,2))
# plt.plot(s15_K60_R1000['mass'] / m_sol, csm_rho)
# plt.ylabel('log rho')
# plt.xlabel('mass (in M_sol)')
# plt.title('CSM')
# plt.yscale('log')

plt.figure(figsize=(6,4))
plt.plot(s15['mass'] / m_sol, np.log10(s15['rho']), marker='o', label='15 $M_{\odot}$', zorder=2, alpha=0.1)
plt.plot(s15_K60_R1000['mass'] / m_sol, np.log10(s15_K60_R1000['rho']), marker='o',
         label='15 $M_{\odot}$ + CSM', zorder=1, alpha=0.8)
plt.plot(s15_K60_R1000['mass'] / m_sol, np.log10(csm_rho), color='black', label='CSM')
plt.ylabel('Log  \u03C1 (gram/cm)')
plt.xlabel('$M_{\odot}$')
plt.title('Mass density profile with vs. without CSM, over mass coordinates')
plt.legend()
plt.savefig('s15_K60_R1000_rho_mass_profile.png')
plt.savefig('s15_K60_R1000_rho_mass_profile.svg')

plt.figure(figsize=(6,4))
plt.plot(s15['radius'] / r_sol, np.log10(s15['rho']), marker='o', label='15 $M_{\odot}$', zorder=2, alpha=0.1)
plt.plot(s15_K60_R1000['radius'] / r_sol, np.log10(s15_K60_R1000['rho']), marker='o',
         label='15 $M_{\odot}$ + CSM', zorder=1, alpha=0.8)
plt.plot(s15_K60_R1000['radius'] / r_sol, np.log10(csm_rho), color='black', label='CSM')
plt.ylabel('Log  \u03C1 (gram/cm)')
plt.xlabel('$R_{\odot}$')
plt.title('Mass density profile with vs. without CSM, over radius coordinates')
plt.legend()
plt.savefig('s15_K60_R1000_rho_radius_profile.png')
plt.savefig('s15_K60_R1000_rho_radius_profile.svg')

plt.show()





