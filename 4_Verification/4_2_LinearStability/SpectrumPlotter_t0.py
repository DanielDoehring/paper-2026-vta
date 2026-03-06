import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

Reals = []
Imags = []


#EigValFileName = "./Eigenvalues_t0_AVI.txt"
#EigValFileName = "./Eigenvalues_t0_WF.txt"
EigValFileName = "./Eigenvalues_t0_FD.txt"

with open(EigValFileName, 'r') as EigValFile:
  Lines = EigValFile.readlines()
  for Line in Lines:
    Entries = Line.split()
    if len(Entries) > 0:  # Exclude blank lines
      imag_part = float(Entries[2][:-1]) # Remove 'i' character at end
      if imag_part >= 0:  # Only include eigenvalues with non-negative imaginary part
        Reals.append(float(Entries[0]))
        Imags.append(imag_part)   

# Plot Section

RWTH_Blue_RGB   = [(0, 84/256, 159/256)]
RWTH_Orange_RGB = [(246/256, 169/256, 0)]
RWTH_Green_RGB  = [(70/256, 171/256, 39/256)]
RWTH_Red_RGB    = [(204/256, 7/256, 30/256)]
RWTH_Petrol_RGB = [(0/256, 152/256, 161/256)]
RWTH_Yellow_RGB = [(225/256, 237/256, 0/256)]
RWTH_Purple_RGB = [(97/256, 33/256, 88/256)]
RWTH_Carmine_RGB = [(161/256, 16/256, 53/256)]

GoldenRatio = (1 + 5 ** 0.5) / 2

InchesX = 4.33071 # = 11 cm = 0.75 * textwidth (in current document)
#InchesX = 5 # I like this
#InchesY = InchesX / GoldenRatio
InchesY = InchesX
#fig, ax = plt.subplots(figsize = (InchesX, InchesY) )
fig, ax = plt.subplots()

ax.set_xlabel("Re")
ax.set_ylabel("Im")

#plt.scatter(Reals, Imags, color = RWTH_Blue_RGB, s = 8.0)
#plt.scatter(Reals, Imags, color = RWTH_Orange_RGB, s = 8.0)
plt.scatter(Reals, Imags, color = RWTH_Carmine_RGB, s = 8.0)

print("Max Real Part: ", max(Reals))

ax.grid(axis ='y', which='both', alpha=0.1, linewidth = 1.5, color ='black')
# Make zero line better visible
gridlines = ax.yaxis.get_gridlines()
#gridlines[0].set_alpha(1)

#ax.set_yticks([0, 0.5, 1.0])

ax.grid(axis ='x', which='both', alpha=0.1, linewidth = 1.5, color ='black')
gridlines = ax.xaxis.get_gridlines()

# Bold imaginary axis
gridlines[4].set_alpha(1)

gridlines = ax.yaxis.get_gridlines()

# Bold real axis
gridlines[0].set_alpha(1)

ax.set_axisbelow(True)  # Hide grid 

#plt.legend(loc='upper left', framealpha=0.7)

#plt.title(r"${\bf \sigma}(J)$ of Adaptive VT/VI")
#plt.title(r"${\bf \sigma}(J)$ of WF VT/VI")
plt.title(r"${\bf \sigma}(J)$ of FD VT/VI")

ax.set_xlim([-4, 4])
ax.set_ylim([0, 140])

#'''
# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)
#'''

plt.tight_layout()  # Title, labels, ... to screen

#plt.savefig("Spectrum_t0_AV.pgf", bbox_inches = 'tight', pad_inches = 0)
#plt.savefig("Spectrum_t0_WF.pgf", bbox_inches = 'tight', pad_inches = 0)
plt.savefig("Spectrum_t0_FD.pgf", bbox_inches = 'tight', pad_inches = 0)

plt.show()
