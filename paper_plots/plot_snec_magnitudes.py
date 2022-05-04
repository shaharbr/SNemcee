import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

filepath = os.path.join('..', 'all_snec_data', 'M15.0_Ni0.02_E0.5_Mix2.0_R0_K0', 'magnitudes.dat')

mag = pd.read_csv(filepath, names=['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'], sep=r'\s+')
# convert seconds to days
mag['time'] = mag['time'] / 86400

fig, ax = plt.subplots()
for filt in ['PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
    ax.plot(mag['time'], mag[filt], label=filt)
    ax.invert_yaxis()
    ax.legend()
    ax.set_title('M15.0_Ni0.02_E0.5_Mix2.0_R0_K0')


plt.show()

print(mag)