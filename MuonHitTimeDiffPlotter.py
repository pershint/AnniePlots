import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as scp
import numpy as np
import sys

NBINS = 150.0

f = uproot.open("./PromptTrig_All_Comb.root")
#f = uproot.open("./MCT_Qgt5_Comb.root")
ftree = f.get("phaseII")
ftree.items()
digitT = ftree.get("digitT")
TrueVtxTime = ftree.get("trueVtxTime")
TrueVtxX = ftree.get("trueVtxX")
TrueVtxY = ftree.get("trueVtxY")
TrueVtxZ = ftree.get("trueVtxZ")
digitX = ftree.get("digitX")
digitY = ftree.get("digitY")
digitType = ftree.get("digitType")
digitZ = ftree.get("digitZ")
evn = ftree.get("eventNumber")
evnums = evn.array()
diT = digitT.array()
diX = digitX.array()
diY = digitY.array()
diZ = digitZ.array()
truX = TrueVtxX.array()
truY = TrueVtxY.array()
truZ = TrueVtxZ.array()
trueT = TrueVtxTime.array()
print("ARRAY OF TRUE MUON VTX TIMES: " + str(trueT))
diType = digitType.array()
tot_entries = len(evn) 
PMT_SOL_residuals = np.array([])
PMT_SOL_wat_residuals = np.array([])
LAPPD_SOL_residuals = np.array([])
LAPPD_SOL_wat_residuals = np.array([])
for e in xrange(tot_entries):
    #Want to calculate all hit residuals based on if
    #Light is moving at C or at C/1.33
    typedata = diType[e]
    pmtind = np.where(typedata==0)[0]
    lappdind = np.where(typedata==1)[0]
    this_TT = trueT[e]
    this_diX = np.array(diX[e])
    this_diY = np.array(diY[e])
    this_diZ = np.array(diZ[e])
    this_diT = np.array(diT[e])
    this_TX = np.array(truX[e])
    this_TY = np.array(truY[e])
    this_TZ = np.array(truZ[e])
    #Calculate speed of light residuals
    L = np.sqrt((this_diX-this_TX)**2 + (this_diY - this_TY)**2 + (this_diZ - this_TZ)**2)
    #Residuals assuming speed of light in vacuum (cm/ns)
    res_SOL = this_diT - (L*1.0/29.97) - this_TT
    pmt_res_SOL = res_SOL[pmtind]
    PMT_SOL_residuals = np.append(pmt_res_SOL,PMT_SOL_residuals)
    lappd_res_SOL = res_SOL[lappdind]
    LAPPD_SOL_residuals = np.append(lappd_res_SOL,LAPPD_SOL_residuals)
    
    #Residuals assuming speed of light in water (cm/ns)
    res_SOL_wat = this_diT - (L*1.33/29.97) - this_TT
    pmt_res_SOL_wat = res_SOL_wat[pmtind]
    PMT_SOL_wat_residuals = np.append(pmt_res_SOL_wat,PMT_SOL_wat_residuals)
    lappd_res_SOL_wat = res_SOL_wat[lappdind]
    LAPPD_SOL_wat_residuals = np.append(lappd_res_SOL_wat,LAPPD_SOL_wat_residuals)

#First, we make a dataset with all PMT-muon times,
#And one with all LAPPD-muon hit times
LAPPDt = np.array([])
PMTt = np.array([])
LAPPDdiff = np.array([])
PMTdiff = np.array([])
mut = np.array([])

for i in range(0,len(evnums)):
    muon_time = trueT[i]
    hit_times = np.array(diT[i])
    digits = diType[i]
    PMT_ind = np.where(digits==0)[0]
    LAPPD_ind = np.where(digits==1)[0]
    pmt_ts = hit_times[PMT_ind]
    lappd_ts = hit_times[LAPPD_ind]
    pmt_diffs = hit_times[PMT_ind] - muon_time
    lappd_diffs = hit_times[LAPPD_ind] - muon_time
    LAPPDdiff=np.append(LAPPDdiff,lappd_diffs)
    PMTdiff=np.append(PMTdiff,pmt_diffs)
    LAPPDt=np.append(LAPPDt,lappd_ts)
    PMTt=np.append(PMTt,pmt_ts)
    mut = np.append(mut,muon_time)

mutimebinedges = np.arange(-5.0, 20.0, 25.0/NBINS)
ltimebinedges = np.arange(-5.0, 30.0, 35.0/NBINS)
ptimebinedges = np.arange(-10.0, 50.0, 60.0/NBINS)
lbinedges = np.arange(-5.0, 25.0, 30.0/NBINS)
pbinedges = np.arange(-40.0, 40.0, 80.0/NBINS)
#Make this into a numpy histogram and we do a fit to the histogram itself
pmt_timehist, pmt_binedges = np.histogram(PMTdiff,bins=pbinedges)
lappd_timehist, lappd_binedges = np.histogram(LAPPDdiff,bins=lbinedges)

sns.set_style("whitegrid")
sns.axes_style("darkgrid")

xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(LAPPDt,bins=ltimebinedges)
plt.ylabel("Entries")
plt.xlabel("LAPPD hit times (ns)")
plt.title("LAPPD hit times after being loaded with DigitBuilder and \n" +\
        "saved to PhaseIITreeMaker",fontsize=30)
plt.show()

xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(PMTt,bins=ptimebinedges)
plt.ylabel("Entries")
plt.xlabel("PMT hit times (ns)")
plt.title("PMT hit times after being loaded with DigitBuilder and \n" +\
        "saved to PhaseIITreeMaker",fontsize=30)
plt.show()

xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(mut,bins=mutimebinedges)
plt.ylabel("Entries")
plt.xlabel("Muon start times (ns)")
plt.title("Muon start times after being loaded with DigitBuilder and \n" +\
        "saved to PhaseIITreeMaker",fontsize=30)
plt.show()

xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(LAPPDdiff,bins=lappd_binedges)
plt.ylabel("Entries")
plt.xlabel("LAPPD hit time - muon time (ns)")
plt.title("Difference in time between LAPPD hits and muon vertex time \n" +\
        "after being loaded into ToolAnalysis",fontsize=30)
plt.show()

xkcd_colors = ['slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(PMTdiff,bins=pmt_binedges)
plt.ylabel("Entries")
plt.xlabel("PMT hit time - muon time (ns)")
plt.title("Difference in time between PMT hits and muon vertex time \n" +\
        "after being loaded into ToolAnalysis",fontsize=30)
plt.show()

xkcd_colors = ['slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#plt.hist(PMT_SOL_wat_residuals,bins=pmt_binedges,label="n=1.33")
plt.hist(PMT_SOL_residuals,bins=pmt_binedges,label="n=1.0")
plt.legend()
plt.ylabel("Entries")
plt.xlabel("PMT hit time residuals (ns)")
plt.title("Hit residual assuming point vertex with different n \n" +\
        "after being loaded into ToolAnalysis",fontsize=30)
plt.show()

xkcd_colors = ['slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#plt.hist(LAPPD_SOL_wat_residuals,bins=lappd_binedges,label="n=1.33")
plt.hist(LAPPD_SOL_residuals,bins=lappd_binedges,label="n=1.0")
plt.legend()
plt.ylabel("Entries")
plt.xlabel("LAPPD hit time residuals (ns)")
plt.title("Hit residual assuming point vertex with different n \n" +\
        "after being loaded into ToolAnalysis",fontsize=30)
plt.show()
