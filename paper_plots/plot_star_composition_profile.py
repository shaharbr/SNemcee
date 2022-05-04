import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

Mzams = 15.0
Ni_mass = 0.07
Ni_mix = 3.0

first_line = ['ncells','nisotopes']
second_line = ['1.0d0 1.0d0 4.0d0 12.0d0 16.0d0 20.0d0 24.0d0 28.0d0 32.0d0 36.0d0 40.0d0 44.0d0 48.0d0 52.0d0 56.0d0 ']
second_line_elements = ['H_p', 'H', 'He', 'C', 'O', 'Ne', 'Mg', 'Si', 'S', 'Cl', 'Ca', 'Sc', 'Ti', 'Cr', 'Fe']
third_line = ["0.0d0 1.0d0 2.0d0 6.0d0 8.0d0 10.0d0 12.0d0 14.0d0 16.0d0 18.0d0 20.0d0 22.0d0 24.0d0 26.0d0 28.0d0 "]
fourth_line = ["cell mass","cell radius"]+["fraction_"+second_line_elements[i] for i in range(len(second_line_elements))]

s17K40_path = os.path.join('data', 'example_star_profiles', 's17.0_K40_R3000.iso.dat')
s17K40 = pd.read_csv(s17K40_path, skiprows=3, sep=r'\s+', names=fourth_line)

s17_path = os.path.join('data', 'example_star_profiles', 's17.0.iso.dat')
s17 = pd.read_csv(s17_path, skiprows=3, sep=r'\s+', names=fourth_line)


R_sol = 6.957* (10**10) # in cm
M_sol = 1.98847 * (10**33)  # in gram
Ni_frac = Ni_mass / Ni_mix

def plot_composition(model, title):
    plt.figure(figsize=(8,6))
    plt.xticks([0, 5, 10, 15])
    plt.minorticks_on()
    plt.xlim(0, 19)
    plt.grid(which='both', alpha=0.2)

    plt.plot(model['cell mass'] / M_sol, model['fraction_H'], label='H', color='red')
    plt.plot(model['cell mass'] / M_sol, model['fraction_He'], label='He', color='blue')
    plt.plot(model['cell mass'] / M_sol, model['fraction_C'], label='C', color='purple')
    plt.plot(model['cell mass'] / M_sol, model['fraction_O'], label='O', color='pink')
    plt.plot([0, Ni_mix, Ni_mix+0.0001, Mzams], [20 * Ni_frac, 20 * Ni_frac, 0, 0], label='Ni x 20', color='cyan')
    plt.title(title)

    plt.ylabel('mass fraction')
    plt.xlabel('Mass (in M_sol)')
    # plt.xscale('log')
    # plt.yscale('log')

model_path = os.path.join('data', 'example_star_profiles',
                                  's'+str(Mzams)+'.iso.dat')
model = pd.read_csv(model_path, skiprows=3, sep=r'\s+', names=fourth_line)
star_max_M = np.max(model['cell mass'] / M_sol)

for K in [60]:
    for R in [1000]:
        model_path = os.path.join('data', 'example_star_profiles',
                                  's'+str(Mzams)+'_K'+str(K)+'_R'+str(R)+'.iso.dat')
        model = pd.read_csv(model_path, skiprows=3, sep=r'\s+', names=fourth_line)
        plot_composition(model, 'Star composition with vs without CSM')
        plt.axvline(star_max_M, linestyle='--', color='black', alpha=0.5, label='original 15 $M_{\odot}$ star')
        plt.xlim(0, star_max_M+3)
        plt.legend(loc='upper right')
        plt.savefig('composition_profile_s'+str(Mzams)+'_K'+str(K)+'_R'+str(R)+'.png')
        plt.savefig('composition_profile_s' + str(Mzams) + '_K' + str(K) + '_R' + str(R) + '.svg')


# model_path = os.path.join('data', 'example_star_profiles',
#                                   's'+str(Mzams)+'.iso.dat')
# model = pd.read_csv(model_path, skiprows=3, sep=r'\s+', names=fourth_line)
# plot_composition(model, str(Mzams)+'_K0_R0')
# plt.savefig('composition_profile_s'+str(Mzams)+'_K0_R0.png')
# plt.savefig('composition_profile_s' + str(Mzams) +'_K0_R0.svg')


plt.show()