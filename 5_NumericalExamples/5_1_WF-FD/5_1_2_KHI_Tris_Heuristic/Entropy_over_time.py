import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Read the data, skipping the first row (header)
df_v = pd.read_csv('./analysis_AVI.dat', delim_whitespace=True, skiprows=1, header=None)
df_FD = pd.read_csv('./analysis_FD.dat', delim_whitespace=True, skiprows=1, header=None)
df_WF = pd.read_csv('./analysis_WF.dat', delim_whitespace=True, skiprows=1, header=None)

RWTH_Blue_RGB   = [(0, 84/256, 159/256)]
RWTH_Orange_RGB = [(246/256, 169/256, 0)]
RWTH_Green_RGB  = [(70/256, 171/256, 39/256)]
RWTH_Red_RGB    = [(204/256, 7/256, 30/256)]
RWTH_Petrol_RGB = [(0/256, 152/256, 161/256)]
RWTH_Yellow_RGB = [(225/256, 237/256, 0/256)]
RWTH_Purple_RGB = [(97/256, 33/256, 88/256)]
RWTH_Carmine_RGB = [(161/256, 16/256, 53/256)]

InchesX = 4.33071 # = 11 cm
'''
GoldenRatio = (1 + 5 ** 0.5) / 2
InchesY = InchesX / GoldenRatio
fig, ax = plt.subplots(figsize = (InchesX, InchesY) )
'''
fig, ax = plt.subplots()

# Plot standard case for each p
plt.plot(df_v[1], df_v.iloc[:, -1] - df_v.iloc[0, -1], color = RWTH_Blue_RGB[0], label = r"Adaptive WF--FD")
#plt.scatter(df_p2[1], df_p2.iloc[:, -1], color = RWTH_Blue_RGB[0])

plt.plot(df_FD[1], df_FD.iloc[:, -1] - df_FD.iloc[0, -1], color = RWTH_Carmine_RGB[0], label = r"Flux--Diff.")
#plt.scatter(df_p3[1], df_p3.iloc[:, -1], color = RWTH_Green_RGB[0])

plt.plot(df_WF[1], df_WF.iloc[:, -1] - df_WF.iloc[0, -1], color = RWTH_Orange_RGB[0], label = r"Weak Form")
#plt.scatter(df_p4[1], df_p4.iloc[:, -1], color = RWTH_Orange_RGB[0])

#plt.title(r"$\Delta S(t)$")
plt.title(r"$S(t) - S(t_0)$")


# Turn on scientific notation
#ax.set_yticks([0, -0.01, -0.02, -0.03, -0.04])

plt.legend(loc = 'lower left')

ax.grid(axis ='y', which='both', alpha=0.1, linewidth = 1.5, color ='black')
ax.grid(axis ='x', which='both', alpha=0.1, linewidth = 1.5, color ='black')
ax.set_axisbelow(True)  # Hide grid 

eps_x = 0.1
ax.set_xlim([-eps_x, 4.6 + eps_x])
#ax.set_ylim([-0.06, 0.002])

ax.set_xlabel(r'$t$')

width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen

plt.savefig('Entropy_KHI_Tris_Heuristic.pgf', bbox_inches = 'tight', pad_inches = 0)

plt.show()
