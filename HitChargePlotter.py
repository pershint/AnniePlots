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

f = uproot.open("./MCT_Qgt5_Comb.root")
ftree = f.get("phaseII")
ftree.items()
digitQ = ftree.get("digitQ")
TrueVtxTime = ftree.get("trueVtxTime")
digitType = ftree.get("digitType")
evn = ftree.get("eventNumber")
evnums = evn.array()
diQ = digitQ.array()
diTypes = digitType.array()
EventHitQs = diQ[ENTRYNUM]
EventdiType = diTypes[ENTRYNUM]
print("FIRSTTIMES: " + str(EventHitQs))
print("DIGIT TYPES: " + str(EventdiType))
print(len(EventHitQs))
evnum = evnums[ENTRYNUM]
binwidth = float(max(diQ[ENTRYNUM]) - min(diQ[ENTRYNUM]))/NBINS

PMT_ind = np.where(EventdiType==0)[0]
LAPPD_ind = np.where(EventdiType==1)[0]

PMTQ = EventHitQs[PMT_ind]
LAPPDQ = EventHitQs[LAPPD_ind]

binedges = np.arange(-5.0, 15.0, 20.0/NBINS)
sns.set_style("whitegrid")
sns.axes_style("darkgrid")

xkcd_colors = [ 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(LAPPDQ,bins=binedges, label="LAPPD Hit Charges")
plt.hist(PMTQ,bins=binedges, label="PMT Hit Charges")
plt.ylabel("Entries")
plt.xlabel("Hit charge")
leg = ax.legend(loc=1)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title(r"Hit charge distribution for event number %i"%(evnum))
plt.show()
