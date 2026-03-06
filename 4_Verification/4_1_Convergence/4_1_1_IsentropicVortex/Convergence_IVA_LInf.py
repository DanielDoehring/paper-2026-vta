#!/usr/bin/python3

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

N = [16.0, 32.0, 64.0, 128.0, 256.0]
LInfErr_AVI = [9.67e-03, 9.58e-04, 1.94e-04, 1.84e-05, 1.82e-06]

LInfErr_WF = [1.78e-02, 2.43e-03, 3.84e-04, 2.49e-05, 2.06e-06]

LInfErr_FD = [6.18e-03, 7.58e-04, 1.18e-04, 8.10e-06, 8.19e-07]

GoldenRatio = (1 + 5 ** 0.5) / 2

InchesX = 4.33071 # = 11 cm = 0.75 * textwidth (in current document)
InchesY = InchesX / GoldenRatio
#fig, ax = plt.subplots(figsize = (InchesX, InchesY) )
fig, ax = plt.subplots()

RWTH_Blue_RGB   = [(0, 84/256, 159/256)]
RWTH_Orange_RGB = [(246/256, 169/256, 0)]
RWTH_Green_RGB  = [(70/256, 171/256, 39/256)]
RWTH_Red_RGB    = [(204/256, 7/256, 30/256)]
RWTH_Petrol_RGB = [(0/256, 156/256, 161/256)]
RWTH_Yellow_RGB = [(225/256, 237/256, 0/256)]
RWTH_Purple_RGB = [(97/256, 33/256, 88/256)]
RWTH_Carmine_RGB = [(161/256, 16/256, 53/256)]

### ACTUAL PLOTTING ###

ax.scatter(N, LInfErr_AVI, label = r'Adaptive', color = RWTH_Blue_RGB, zorder = 3)
ax.plot(N, LInfErr_AVI, color = RWTH_Blue_RGB[0], linestyle='dashed', zorder = 3)

ax.scatter(N, LInfErr_WF, label = r'Weak Form', color = RWTH_Orange_RGB, zorder = 2)
ax.plot(N, LInfErr_WF, color = RWTH_Orange_RGB[0], linestyle='dashed', zorder = 2)

ax.scatter(N, LInfErr_FD, label = r'Flux--Diff.', color = RWTH_Carmine_RGB, zorder = 2)
ax.plot(N, LInfErr_FD, color = RWTH_Carmine_RGB[0], linestyle='dashed', zorder = 2)

ax.loglog(N, np.multiply(6e3, np.power(N, -4) ), linestyle='dashed', zorder = 1,
            label = r'$\mathcal{O}\left(\Delta x^4\right)$',
            color = 'black') # Order two line fitted

# Turn on logscale (no native support for logarithmic scatter)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$N_x = N_y$')

### LIMITS SECTION ###
eps_x = 0.1
ax.set_xlim([15, 280])

ax.set_ylim([5e-7, 0.03])

### GRID SECTION ###
ax.grid(axis ='both', which='major', alpha=0.1, linewidth = 1.5, color ='black')
ax.set_axisbelow(True)  # Hide grid behind bars

### LEGEND SECTION ###
ax.legend(loc = "lower left")

### TICKS SECTION ###

ax.set_xticks(N)

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.get_xaxis().set_tick_params(which='minor', size=0)
ax.get_xaxis().set_tick_params(which='minor', width=0) 

ax.set_xticklabels([r"$16$", r"$32$", r"$64$", r"$128$", r"$256$"])

# Make bounding lines thicker
'''
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
'''

### TITLE SECTION ###
plt.title(r"$e^{L^\infty}_\rho$ (Isentropic Vortex)")

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen
plt.savefig('Convergence_IVA_LInf.pgf', bbox_inches = 'tight', pad_inches = 0)
plt.show()
