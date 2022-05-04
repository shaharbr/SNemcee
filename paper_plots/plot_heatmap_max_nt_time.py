import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

cd .
df_last_time = pd.read_csv('../last_nt_time.csv', usecols=['Mzams', 'E_final', 'last_time'])
df_last_time['last_time'] = df_last_time['last_time'] / 86400
df_last_time = df_last_time.pivot('Mzams', 'E_final', 'last_time')

df_last_nt = pd.read_csv('../last_nt_time.csv', usecols=['Mzams', 'E_final', 'last_nt'])
df_last_nt = df_last_nt.pivot('Mzams', 'E_final', 'last_nt')

def mix_palette():
    palette = sns.color_palette("viridis", 30)
    palette[0] = 'lightgray'
    return palette

# TODO fill up the non-calculated boxes with predicted values (4M etc), so that the pixel sizes will be proportional to size.
#  also take note that the below-border pixels should be different colors - lighter purple etc (smear based on closes neighbor)
#  then add the labled on the axes to only show the values actually calculated (so we're not lying)
# for Mzams in [9.5, 10.5, 12.5, 13.5, 14., 14.5, 15.5, 16., 16.5]:
# for E in [0.45, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.95, 1.,
#           1.05, 1.1, 1.15, 1.2, 1.25, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65]:


plt.figure()
ax1 = sns.heatmap(df_last_time)
plt.suptitle('Days calculated (Max number of steps allowed: 4 million)', fontsize=16)
plt.figure()
ax2 = sns.heatmap(df_last_nt, cmap=mix_palette(), square=True,
                  xticklabels = True, yticklabels = True, norm=LogNorm())

plt.suptitle('Steps calculated until day 200 or max steps reached (4 million)', fontsize=16)

plt.show()
