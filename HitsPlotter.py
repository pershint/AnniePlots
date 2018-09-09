import uproot
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('poster')

f = uproot.open("./RecoTree.root")
ftree = f.get("phaseII")
ftree.items()
fTVZ = ftree.get("trueVtxZ")
fTVY = ftree.get("trueVtxY")
fTVX = ftree.get("trueVtxX")
trueXvtx = fTVX.array()
trueZvtx = fTVZ.array()
trueYvtx = fTVY.array()
digitX = ftree.get("digitX")
diX = digitX.array()
digitY = ftree.get("digitY")
diY = digitY.array()
digitZ = ftree.get("digitZ")
diZ = digitZ.array()
sns.set_style("whitegrid")
sns.axes_style("darkgrid")
xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(diZ[0],diX[0], linestyle='none', marker='o',
markersize=4, label="Digit Positions")
ax.plot(trueZvtx[0], trueXvtx[0], linestyle='none', 
markersize=20, label="True Vertex", marker='*')
plt.ylabel("X [cm]")
plt.xlabel("Z [cm]")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title(r"Hits in ANNIE from simulated $\nu_{\mu}$ interaction producing a muon")
plt.show()
