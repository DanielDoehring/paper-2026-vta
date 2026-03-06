import matplotlib.pyplot as plt
import math # for pi

TimeStampsAVI = []
EnstropyAVI = []

DataFileName = "./analysis_AVI_k3_32.dat"
with open(DataFileName, 'r') as DataFile:
  Lines = DataFile.readlines()
  isFirstLine = True
  for Line in Lines:
    if isFirstLine:
      isFirstLine = False
      continue
    Entries = Line.split()
    if len(Entries) > 0:  # Exclude blank lines
      TimeStampsAVI.append(float(Entries[1]))
      EnstropyAVI.append(float(Entries[6]))

TimeStampsEC = []
EnstropyEC = []

DataFileName = "./analysis_EC_k3_32.dat"
with open(DataFileName, 'r') as DataFile:
  Lines = DataFile.readlines()
  isFirstLine = True
  for Line in Lines:
    if isFirstLine:
      isFirstLine = False
      continue
    Entries = Line.split()
    if len(Entries) > 0:  # Exclude blank lines
      TimeStampsEC.append(float(Entries[1]))
      EnstropyEC.append(float(Entries[6]))

TimeStampsWang512 = []
EnstrophyWang512 = []

mu = 1/1600

DataFileName = "./Wang_512.txt"
with open(DataFileName, 'r') as DataFile:
  Lines = DataFile.readlines()
  for Line in Lines:
    Entries = Line.split()
    if len(Entries) > 0:  # Exclude blank lines
      TimeStampsWang512.append(float(Entries[0]))
      EnstrophyWang512.append(float(Entries[1]) / (2*mu))

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
InchesY = InchesX / GoldenRatio
#fig, ax = plt.subplots(figsize = (InchesX, InchesY) )
fig, ax = plt.subplots()

ax.set_xlabel(r"$t$")
#ax.set_ylabel(r"$\bar{\epsilon}$ ", rotation = 0) # Space behind rho such that it does not overlap/is too close with y-axis
ax.set_ylabel(r"$\bar{\epsilon}$ ")


ax.plot(TimeStampsWang512, EnstrophyWang512, color = RWTH_Green_RGB[0], label = "Reference")

ax.plot(TimeStampsAVI, EnstropyAVI, color = RWTH_Blue_RGB[0], label = "Adaptive")
ax.plot(TimeStampsEC, EnstropyEC, color = RWTH_Carmine_RGB[0], label = "Flux--Diff.")

ax.grid(axis ='y', which='both', alpha=0.1, linewidth = 1.5, color ='black')
ax.grid(axis ='x', which='both', alpha=0.1, linewidth = 1.5, color ='black')

ax.set_axisbelow(True) # Hide grid behind data

ax.set_xlim([-0.5, 20.5])

plt.legend(loc='upper right', framealpha=0.9)
plt.title("Enstrophy")

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen

Ending = ".pgf"
plt.savefig("Enstrophy" + Ending, bbox_inches = 'tight', pad_inches = 0)

plt.show()
