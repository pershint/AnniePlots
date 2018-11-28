import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as scp
import numpy as np
import sys

sns.set_context('poster')
try:
    ENTRYNUM = int(sys.argv[1])
except:
    print("Something went wrong with your given event to plot")
    print("Setting entry number to plot as 0")
    ENTRYNUM = 0

NBINS = 30
#####/TUNABLES FOR GRAPH OUTPUT####

f = uproot.open("./AllGridEffInd.root")
ftree = f.get("phaseII")
ftree.items()
svFOM = ftree.get("seedVtxFOM")
TrueVtxTime = ftree.get("trueVtxTime")
digitType = ftree.get("digitType")
evn = ftree.get("eventNumber")
evnums = evn.array()
theFOMs = svFOM.array()
diTypes = digitType.array()
EventFOMs = theFOMs[ENTRYNUM]
EventdiType = diTypes[ENTRYNUM]
print("FIRSTTIMES: " + str(EventFOMs))
print("DIGIT TYPES: " + str(EventdiType))
print(len(EventFOMs))
evnum = evnums[ENTRYNUM]
binwidth = float(max(theFOMs[ENTRYNUM]) - min(theFOMs[ENTRYNUM]))/NBINS

binedges = np.arange(0, 100.0, 100.0/NBINS)
sns.set_style("whitegrid")
sns.axes_style("darkgrid")

xkcd_colors = [ 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(EventFOMs,bins=binedges)
plt.ylabel("Entries")
plt.xlabel("Point Position FOMs")
plt.title(r"Grid Seed Vtx FOM distribution for event number %i"%(evnum))
plt.show()
