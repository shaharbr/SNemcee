import os
import math
import pandas as pd
import numpy as np


'''
This code was adapted from the original wind adding code written by K. Azalee Bostroem : github.com/abostroem/SNEC
'''

def add_wind(K, R, Mzams, dir):
    KK = K *(10**17)
    msol = 1.98e33
    rsol = 6.96e10
    ggrav = 6.6742e-8
    velocity_of_wind = 1.0e6
    delta_r = 1  # Radius step in solar radii

    ### --------------------- Parameters --------------------

    mainfolder = dir+'profiles/'


    rho_attach_gl = []
    vel_esc_gl = []
    mdot_gl = []
    mass_in_wind_gl = []
    tau_wind_gl = []

    fname = mainfolder+'s{}.short'.format(str(Mzams))  # change
    fname_iso = mainfolder + 's{}.iso.dat'.format(str(Mzams))  # change

    profile_df = pd.read_csv(fname, skiprows=1, names=['index', 'mass', 'radius', 'temp', 'rho', 'vel'], sep='\s+')
    tot_rad = np.max(profile_df['radius']) / rsol
    added_R = R
    R += tot_rad


    ### --------------------- Read the profile of the core --------------------

    mass = []
    radius = []
    temp = []
    rho = []
    rho_log = []
    vel = []

    for l in open(fname, 'r').readlines():
        s = l.split()
        if len(s) != 1:
            mass.append(float(s[1]) / msol)
            radius.append(float(s[2]) / rsol)
            temp.append(float(s[3]))
            rho.append(float(s[4]))
            rho_log.append(math.log10(float(s[4])))
            vel.append(float(s[5]))
    # Attach to profile when densities are the same
    rho_wind_attach = float(KK) / ((radius[-1] * rsol) ** 2)

    vel_esc = math.sqrt(2 * ggrav * mass[-1] * msol / (radius[-1] * rsol))
    mdot_msol_yr_max = rho_wind_attach * 4 * math.pi * (radius[
                                                            -1] * rsol) ** 2 * velocity_of_wind * 365 * 24 * 3600 / msol

    rho_attach_gl.append(rho_wind_attach)
    vel_esc_gl.append(vel_esc)
    mdot_gl.append(mdot_msol_yr_max)

    ### -------------------- Building profiles with the wind ----------------------

    radius_wind = []  # in solar radii
    rho_wind = []
    rho_log_wind = []
    mass_wind = []  # in solar masses
    temp_wind = []
    vel_wind = []

    for l in range(10000000):  # (just big number to cover all the profile)
        # if the radius coordinate (=line number in the profile) is within the radius of the star,
        # just copy the content of the star
        if l < len(radius) and rho[l] > rho_wind_attach:
            radius_wind.append(radius[l])
            rho_wind.append(rho[l])
            rho_log_wind.append(rho_log[l])
            mass_wind.append(mass[l])
            temp_wind.append(temp[l])
            vel_wind.append(vel[l])
        # if the radius coordinate is beyond the radius of the star, check if its beyond the chosen radius
        # for the wind, and if it is - break
        else:
            if (radius_wind[l - 1] + delta_r) > float(R):
                break
            # if the radius coordinate is beyond the radius of the star and within the radius of the wind,
            # start calculating and adding there the wind
            radius_wind.append(radius_wind[l - 1] + delta_r)
            rho_wind.append(rho_wind_attach * (radius[-1] / radius_wind[l]) ** 2)
            rho_log_wind.append(math.log10(rho_wind[l]))
            mass_wind.append(mass_wind[l - 1] + 4.0 * math.pi / 3.0 * (radius_wind[l] ** 3 -
                                                                       radius_wind[l - 1] ** 3) *
                             rho_wind[l] * rsol ** 3 / msol)
            temp_wind.append(temp[-1])
            vel_wind.append(velocity_of_wind)

    mass_in_wind = mass_wind[-1] - mass[-1]
    # print('mass in wind', mass_in_wind)
    mass_in_wind_gl.append(mass_in_wind)

    ### ------------ calculating the optical depth of the wind --------------
    # This is informational - SNEC will recalculate more carefully
    kappa = 1.0e-4
    tau_wind = 0
    for l in range(len(radius_wind) - 1):
        if radius_wind[len(radius_wind) - 1 - l - 1] < radius[-1]:
            break
        tau_wind = tau_wind + kappa * rsol * (radius_wind[len(radius_wind) - 1 - l]
                                              - radius_wind[len(radius_wind) - 1 - l - 1]) * rho_wind[
                                  len(radius_wind) - 1 - l - 1]

    tau_wind_gl.append(tau_wind)

    ### -------------------- reading a composition file ------------------------

    iso_lines = []

    n_line = 0

    for l in open(fname_iso, 'r').readlines():
        if n_line == 0:
            s = l.split()
            ncells = int(s[0])
            niso = int(s[1])
        elif n_line == 1:
            A_line = l
        elif n_line == 2:
            Z_line = l
        else:
            iso_lines.append(l)
            s = l.split()
            last_line = s

        n_line = n_line + 1

    ### --------------------- Output the profiles -----------------------------
    # note, that here we convert mass and radius in cgs units, used by the code.
    # in the script they are in the solar masses and radii for convenience of plotting

    outfile = open(os.path.join(mainfolder, 's{}_K{}_R{}.short'.format(str(Mzams), str(K), str(added_R))), "w")

    outfile.write(str(len(radius_wind)) + '\n')

    for l in range(len(radius_wind)):
        if l == len(radius_wind) - 1:
            outfile.write(str(l + 1) + '   ' + str(mass_wind[l] * msol) + '   ' + str(radius_wind[l] * rsol)
                          + '   ' + str(temp_wind[l]) + '   ' + str(rho_wind[l]) + '   '
                          + str(vel_wind[l]))
        else:
            outfile.write(str(l + 1) + '   ' + str(mass_wind[l] * msol) + '   ' + str(radius_wind[l] * rsol)
                          + '   ' + str(temp_wind[l]) + '   ' + str(rho_wind[l]) + '   '
                          + str(vel_wind[l]) + '\n')

    outfile.write('\n')

    outfile.close()
    outfile_iso = open(os.path.join(mainfolder, 's{}_K{}_R{}.iso.dat'.format(str(Mzams), str(K), str(added_R))), "w")

    outfile_iso.write(str(ncells + 1) + '   ' + str(niso) + '\n')
    outfile_iso.write(A_line)
    outfile_iso.write(Z_line)

    for l in range(len(iso_lines)):
        outfile_iso.write(iso_lines[l])
    # outfile_iso.write('\n') # for some reason it was used before...

    outfile_iso.write(str(mass_wind[-1] * msol) + '        ' + str(radius_wind[-1] * rsol) + '          ')

    for l in range(2, len(last_line)):
        outfile_iso.write(last_line[l] + '   ')

    outfile_iso.write('\n')

    outfile_iso.close()

    ######### ------------ table with the wind properties

    outfile_info = open(mainfolder + '/info.dat', "w")

    outfile_info.write(str(Mzams)+ ": rho attach=" + str(rho_attach_gl[0]) + ", vel esc=" +
                           str(vel_esc_gl[0]) + ", mdot=" + str(mdot_gl[0]) + ", mass in wind=" +
                           str(mass_in_wind_gl[0]) + ", tau wind=" + str(tau_wind_gl[0]) + '\n')

    outfile_info.close()
