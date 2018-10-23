import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sns.set_context('poster')

try:
    eventnum = int(sys.argv[1])
except:
    print("Something went wrong with your given event to plot")
    print("Setting entry number to plot as 0")
    eventnum = 0

f = uproot.open("./RecoGridSeed_5LAPPD_Comb.root")
ftree = f.get("phaseII")
ftree.items()
fTVZ = ftree.get("trueVtxZ")
fTVY = ftree.get("trueVtxY")
fTVX = ftree.get("trueVtxX")
fPVZ = ftree.get("pointPosZ")
fPVY = ftree.get("pointPosY")
fPVX = ftree.get("pointPosX")
fRC = ftree.get("recoStatus")
trueXvtx = fTVX.array()
trueZvtx = fTVZ.array()
trueYvtx = fTVY.array()
pointposXvtx = fPVX.array()
pointposZvtx = fPVZ.array()
pointposYvtx = fPVY.array()
recoStatus = fRC.array()
fSVZ = ftree.get("seedVtxZ")
fSVY = ftree.get("seedVtxY")
fSVX = ftree.get("seedVtxX")
seedXvtx = fSVX.array()
seedZvtx = fSVZ.array()
seedYvtx = fSVY.array()
sns.set_style("whitegrid")
sns.axes_style("darkgrid")
xkcd_colors = ['magenta', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))

#Print a warning if recoStatus==0 when all is said and done
if recoStatus[eventnum]!=0:
    print("This event ultimately failed the reconstruction chain somewhere.")
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(seedZvtx[eventnum], seedXvtx[eventnum], linestyle='none', 
markersize=5, label="Vertex Seeds", marker='o')
ax.plot(trueZvtx[eventnum], trueXvtx[eventnum], linestyle='none', 
markersize=20, label="True Vertex", marker='*')
ax.plot(pointposZvtx[eventnum], pointposXvtx[eventnum], linestyle='none', 
markersize=20, label="PointPosition Vertex", marker='*')
plt.ylabel("X [cm]")
plt.xlabel("Z [cm]")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title("Seed generator results from simulated event")
plt.show()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(seedZvtx[eventnum], seedYvtx[eventnum], linestyle='none', 
markersize=5, label="Vertex Seeds", marker='o')
ax.plot(trueZvtx[eventnum], trueYvtx[eventnum], linestyle='none', 
markersize=20, label="True Vertex", marker='*')
ax.plot(pointposZvtx[eventnum], pointposYvtx[eventnum], linestyle='none', 
markersize=20, label="PointPosition Vertex", marker='*')
plt.ylabel("Y [cm]")
plt.xlabel("Z [cm]")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title("Seed generator results from simulated event")
plt.show()
