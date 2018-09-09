import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_context('poster')

f = uproot.open("./RecoTree.root")
ftree = f.get("phaseII")
ftree.items()
fPVZ = ftree.get("pointPosVtxZ")
fPVY = ftree.get("pointPosVtxY")
fPVX = ftree.get("pointPosVtxX")
fPVS = ftree.get("pointPosVtxStatus")
pointPosXvtx = fPVX.array()
print("POINTPOSX: " + str(pointPosXvtx))
pointPosZvtx = fPVZ.array()
pointPosYvtx = fPVY.array()
pointPosSvtx = fPVS.array()
print("POINTPOSSTATUS: " + str(pointPosSvtx))
fTVZ = ftree.get("trueVtxZ")
fTVY = ftree.get("trueVtxY")
fTVX = ftree.get("trueVtxX")
trueXvtx = fTVX.array()
trueZvtx = fTVZ.array()
trueYvtx = fTVY.array()
validfits = np.where(pointPosSvtx==0)[0]
trueXVF = np.array(trueXvtx[validfits])
trueZVF = np.array(trueZvtx[validfits])
pointPosXVF = np.array(pointPosXvtx[validfits])
pointPosZVF = np.array(pointPosZvtx[validfits])

#Plotting
sns.set_style("whitegrid")
sns.axes_style("darkgrid")
xkcd_colors = ['purple', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
print(np.sqrt((pointPosXVF**2) - (trueXVF**2)))
ax.plot(  np.sqrt(((pointPosZvtx-trueZvtx)**2)), np.sqrt(((pointPosXvtx-trueXvtx)**2)), linestyle='none', 
markersize=8, label=r"|true-pointPos|", marker='o')
plt.xlabel(r"|PointPos-True|$_{Z}$")
plt.ylabel(r"|PointPos-True|$_{X}$")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title("Difference in PointPosition fitter and true vertex for XZ coordinates\n"+\
        "Only for events with a generated muon and valid PointPosition")
plt.show()
