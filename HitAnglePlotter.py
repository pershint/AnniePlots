import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as scp
import numpy as np
import sys

NBINS = 150

f = uproot.open("./AllGridEffInd.root")
ftree = f.get("phaseII")
ftree.items()
digitQ = ftree.get("digitQ")
TrueVtxTime = ftree.get("trueVtxTime")
TrueVtxX = ftree.get("trueVtxX")
TrueVtxY = ftree.get("trueVtxY")
TrueVtxZ = ftree.get("trueVtxZ")
TrueDirX = ftree.get("trueDirX")
TrueDirY = ftree.get("trueDirY")
TrueDirZ = ftree.get("trueDirZ")
digitX = ftree.get("digitX")
digitY = ftree.get("digitY")
digitType = ftree.get("digitType")
digitZ = ftree.get("digitZ")
digitT = ftree.get("digitT")
evn = ftree.get("eventNumber")
evnums = evn.array()
diT = digitT.array()
diX = digitX.array()
diY = digitY.array()
diZ = digitZ.array()
diQ = digitQ.array()
truX = TrueVtxX.array()
truY = TrueVtxY.array()
truZ = TrueVtxZ.array()
truDirX = TrueDirX.array()
truDirY = TrueDirY.array()
truDirZ = TrueDirZ.array()
trueT = TrueVtxTime.array()
diType = digitType.array()

#Let's look at the phi distribution for events, because yeah
#Get the phi
allphis_pmt = np.array([])
allQs_pmt = np.array([])
allphis_lappd = np.array([])
allQs_lappd = np.array([])
for e in xrange(len(evn)):
    thisQ = diQ[e]
    typedata = diType[e]
    #For each digit, we need the dir from trueVtx to the digit
    thisDiX =np.array(diX[e])
    thisDiY =np.array(diY[e])
    thisDiZ =np.array(diZ[e])
    thistruX =np.array(truX[e])
    thistruY =np.array(truY[e])
    thistruZ =np.array(truZ[e])
    magdiff = np.sqrt((thisDiX-thistruX)**2 + (thisDiY-thistruY)**2 + (thisDiZ-thistruZ)**2)
    pX = (thisDiX - thistruX)/(magdiff)
    pY = (thisDiY - thistruY)/(magdiff)
    pZ = (thisDiZ - thistruZ)/(magdiff)
    thistrudirX = truDirX[e]
    thistrudirY = truDirY[e]
    thistrudirZ = truDirZ[e]
    phi_rad = np.arccos(pX*thistrudirX + pY*thistrudirY + pZ*thistrudirZ)
    phi_deg = phi_rad * (180/np.pi)
    pmtind = np.where(typedata==0)[0]
    lappdind = np.where(typedata==1)[0]
    allphis_pmt = np.append(phi_deg[pmtind],allphis_pmt)
    allQs_pmt = np.append(thisQ[pmtind],allQs_pmt)
    allphis_lappd = np.append(phi_deg[lappdind],allphis_lappd)
    allQs_lappd = np.append(thisQ[lappdind],allQs_lappd)

sns.set_style("whitegrid")
sns.axes_style("darkgrid")
xkcd_colors =  [ 'warm pink', 'slate blue', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(allphis_lappd,bins=NBINS,label="LAPPDs",normed=True)
plt.hist(allphis_pmt,bins=NBINS,label="PMTs",alpha=0.6,normed=True)
plt.legend(fontsize=20)
for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
plt.ylabel("Probability of hit",fontsize=30)
plt.xlabel("Digit angle from true vertex direction (degrees)",fontsize=30)
plt.title("Anglular distribution of hits from true vertex direction \n"+\
        "GENIE-simulated CCQE events, WCSim detector simulation",fontsize=34)
plt.show()
