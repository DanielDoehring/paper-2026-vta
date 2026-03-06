import matplotlib.pyplot as plt

TimeStampsAVI = []
EKinAVI = []

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
      EKinAVI.append(float(Entries[4]))

TimeStampsEC = []
EKinEC = []

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
      EKinEC.append(float(Entries[4]))

TimeStampsWang512 = []
EKinWang512 = []      

DataFileName = "./EKin_Wang_512.txt"
with open(DataFileName, 'r') as DataFile:
  Lines = DataFile.readlines()
  for Line in Lines:
    Entries = Line.split()
    if len(Entries) > 0:  # Exclude blank lines
      TimeStampsWang512.append(float(Entries[0]))
      EKinWang512.append(float(Entries[1]))

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
ax.set_ylabel(r"$\overline{E}_\mathrm{kin}$")

ax.plot(TimeStampsWang512, EKinWang512, color = RWTH_Green_RGB[0], label = "Reference")
ax.plot(TimeStampsAVI, EKinAVI, color = RWTH_Blue_RGB[0], label = "Adaptive")
ax.plot(TimeStampsEC, EKinEC, color = RWTH_Carmine_RGB[0], label = "Flux--Diff.")

ax.grid(axis ='y', which='both', alpha=0.1, linewidth = 1.5, color ='black')
ax.grid(axis ='x', which='both', alpha=0.1, linewidth = 1.5, color ='black')

ax.set_axisbelow(True) # Hide grid behind data

plt.legend(loc='upper right', framealpha=0.9) # Legend already in enstrophy plot
plt.title("Kinetic Energy")

ax.set_xlim([-0.5, 20.5])

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen

Ending = ".pgf"
plt.savefig("EKin" + Ending, bbox_inches = 'tight', pad_inches = 0)

plt.show()
