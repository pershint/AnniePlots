import uproot
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('poster')
import sys
import numpy as np

def IDsInRange(xdata, ydata, zdata, IDdata, ymin, ymax, thetamin, thetamax):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    theta = []
    for i in xrange(len(xdata)):
        thistheta = None
        x = xdata[i]
        z = zdata[i]
        isnegative = False
        if x <0:
            isnegative = True
            x = -x
        r = np.sqrt(z**2 + x**2)
        if r == 0:
            thistheta = np.arccos(0)*180.0/np.pi
        else:
            thistheta = np.arccos(z/r)*180.0/np.pi 
        if isnegative:
            thistheta = (360.0 - thistheta)
            #thistheta = (180.0 + thistheta)
        #Now, do the transormation to beam theta
        if thistheta < 180:
            thistheta = - thistheta
        else:
            thistheta =  (360.0 - thistheta)

        print(thistheta)
        theta.append(thistheta) #yeah that's confusing labeling
    y = ydata
    idCut = []
    thetaCut = []
    for j,ID in enumerate(IDdata):
        if ydata[j] >ymin and ydata[j] < ymax:
            if theta[j] > thetamin and theta[j] < thetamax:
                idCut.append(ID)
    plt.hist(idCut)
    ax.set_xlabel("LAPPD Theta (deg)", fontsize=22)
    ax.set_ylabel("Y position (m)", fontsize=22)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
    plt.title("LAPPDIDs in y and theta region identified")
    plt.show()

def YVSTheta(xdata, ydata, zdata):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    theta = []
    for i in xrange(len(xdata)):
        thistheta = None
        x = xdata[i]
        z = zdata[i]
        isnegative = False
        if x <0:
            isnegative = True
            x = -x
        r = np.sqrt(z**2 + x**2)
        if r == 0:
            thistheta = np.arccos(0)*180.0/np.pi
        else:
            thistheta = np.arccos(z/r)*180.0/np.pi 
        if isnegative:
            thistheta = (360.0 - thistheta)
            #thistheta = (180.0 + thistheta)
        #Now, do the transormation to beam theta
        if thistheta < 180:
            thistheta = - thistheta
        else:
            thistheta =  (360.0 - thistheta)

        print(thistheta)
        theta.append(thistheta) #yeah that's confusing labeling
    y = ydata
    #X,Y = np.meshgrid(x, y)
    plt.plot(theta,y,linestyle='none',markersize=4, marker='o')
    ax.set_xlabel("LAPPD Theta (deg)", fontsize=22)
    ax.set_ylabel("Y position (m)", fontsize=22)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
    plt.title("Positions of all LAPPDS hit in file loaded")
    plt.show()


f = uproot.open("./RecoGridSeed_5LAPPD_Comb.root")
#GIVE AN EVENT NUMBER TO LOOK AT
try:
    eventnum = int(sys.argv[1])
except:
    print("Something went wrong with your given event to plot")
    print("Setting entry number to plot as 0")
    eventnum = 0
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
digitType = ftree.get("digitType")
diType = digitType.array()
digitID = ftree.get("digitDetID")
diID = digitID.array()
sns.set_style("whitegrid")
sns.axes_style("darkgrid")
xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(diZ[eventnum],diX[eventnum], linestyle='none', marker='o',
markersize=4, label="Digit Positions")
ax.plot(trueZvtx[eventnum], trueXvtx[eventnum], linestyle='none', 
markersize=20, label="True Vertex", marker='*')
plt.ylabel("X [cm]")
plt.xlabel("Z [cm]")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title(r"Hits in ANNIE from simulated $\nu_{\mu}$ interaction producing a muon")
plt.show()

#Plot all hits at front of tank
fronthits = np.where(diZ[eventnum]>0.0)[0]
frontX = diX[eventnum][fronthits]
frontY = diY[eventnum][fronthits]
frontID = diID[eventnum][fronthits]
print("THE FRONT IDS: " + str(frontID))
#Now, only plot hits with a single id
IDcut = 59 
thisid = np.where(frontID==IDcut)[0]
print("THIS ID: " + str(thisid))
frontX = frontX[thisid]
frontY = frontY[thisid]
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(frontX,frontY, linestyle='none', marker='o',
markersize=7, label="Digit Positions")
plt.ylabel("Y [cm]")
plt.xlabel("X [cm]")
leg = ax.legend(loc=2)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title(r"Hits in ANNIE from simulated $\nu_{\mu}$ interaction producing a muon")
plt.show()

YVSTheta(diX[eventnum],diY[eventnum],diZ[eventnum])
IDsInRange(diX[eventnum],diY[eventnum],diZ[eventnum],diID[eventnum],-30,0,-40,-20)
