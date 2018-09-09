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
fSVZ = ftree.get("seedVtxZ")
fSVY = ftree.get("seedVtxY")
fSVX = ftree.get("seedVtxX")
seedXvtx = fSVX.array()
seedZvtx = fSVZ.array()
seedYvtx = fSVY.array()
sns.set_style("whitegrid")
sns.axes_style("darkgrid")
xkcd_colors = ['grass', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(seedZvtx[0], seedXvtx[0], linestyle='none', 
markersize=5, label="Vertex Seeds", marker='o')
ax.plot(trueZvtx[0], trueXvtx[0], linestyle='none', 
markersize=20, label="True Vertex", marker='*')
plt.ylabel("X [cm]")
plt.xlabel("Z [cm]")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title("Seed generator results from simulated event")
plt.show()
