import uproot
from mpl_toolkits.mplot3d import Axes3D
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

f = uproot.open("./AllGridEffInd.root")
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
ax = fig.add_subplot(111, projection='3d')
seedx, seedy, seedz = [], [], []
ax.scatter(seedZvtx[eventnum],seedXvtx[eventnum], seedYvtx[eventnum],label="Vertex seeds")
ax.scatter(trueZvtx[eventnum], trueXvtx[eventnum], trueYvtx[eventnum],label="True Vertex",c='black',s=60,marker='*')
ax.set_xlabel("Z position (mm)", fontsize=22)
ax.set_ylabel("X position (mm)", fontsize=22)
ax.set_zlabel("Y position (mm)", fontsize=22)
for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(20)
for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
plt.title("Seed points colored by their point fit FOM in ANNIE",fontsize=34)
plt.show()

