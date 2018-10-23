import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as scp
import numpy as np
import sys

def Chi2OverNDOF(fit,dat,sig,NumVar):
    #Do not have bins with no data contribute to fit
    dat_nonzero_ind = np.where(dat>0)[0]
    dat = dat[dat_nonzero_ind]
    fit = fit[dat_nonzero_ind]
    #assume uncertainty for each dat bin is sqrt(dat)
    k = float(len(dat))
    c = float(NumVar)+1
    Chi2_eachbin = ((fit-dat)**2)/dat
    return np.sum(Chi2_eachbin)/(k-c)



sns.set_context('poster')
gaussTimesExpo= lambda x,C,l,m,s,t: ((C*np.exp(-(x-t)*l))*((s**2)*2*np.pi)**(-1./2)*np.exp(-1/2*(x-m)**2/s**2))
gaussPlusExpo= lambda x,C,D,l,m,s,t: (C*(np.exp(-(x-t)*l))+D*((s**2)*2*np.pi)**(-1./2)*np.exp(-1/2*(x-m)**2/s**2))
landauPlusGauss= lambda x,C,l1,l2,D,m,s: C*((1./np.sqrt(2*np.pi))*np.exp(-(1./2)*(((x-m)/l1) + np.exp(-(x-m)/l2))))+D*((s**2)*2*np.pi)**(-1./2)*np.exp(-1/2*(x-m)**2/s**2)
#landauPlusExpo = lambda x,H,t,l,C,m,l1,l2: H*(np.exp(-(x-t)*l)) + C*((1./np.sqrt(2*np.pi))*np.exp(-(1./2)*(((x-m)/l1) + np.exp(-(x-m)/l2))))
landau = lambda x,C,m,l1,l2: C*((1./np.sqrt(2*np.pi))*np.exp(-(1./2)*(((x-m)/l1) + np.exp(-(x-m)/l2))))
#wald = lambda x,C,m,s,d:  C * np.sqrt(m / (2.* np.pi *(x**s))) * np.exp(-(m*((x-d)**2))/(2.*(d**2)*x))
#exponential_decay_pcon = lambda x, C, l, a: a + (C*np.exp(-x*l))
    #gauss = lambda x, m, s: (s**2*2*np.pi)**(-1/2)*np.exp(-1/2*(x-m)**2/s**2)

#####TUNABLES FOR GRAPH OUTPUT#####
thefunc = landau
NumVariables = 4
try:
    ENTRYNUM = int(sys.argv[1])
except:
    print("Something went wrong with your given event to plot")
    print("Setting entry number to plot as 0")
    ENTRYNUM = 0

NBINS = 100 
#p0 = [100., 100, 100.]
#p0 = [100.,100.,1.,12.,1.,10] #converges with ENTRYNUM=1 for gaussExpo
#p0 = [10., 100., 1., 10.,1.,1.] #Exploring for landau distribution convergence 
p0 = [100., 10., 1., 1.]
#####/TUNABLES FOR GRAPH OUTPUT####

f = uproot.open("./RecoGridSeed_5LAPPD_Comb.root")
ftree = f.get("phaseII")
ftree.items()
digitT = ftree.get("digitT")
digitType = ftree.get("digitType")
evn = ftree.get("eventNumber")
evnums = evn.array()
diT = digitT.array()
diType = digitType.array()
FirstTimes = diT[ENTRYNUM]
diTypes = diType[ENTRYNUM]
print("FIRSTTIMES: " + str(FirstTimes))
print("DIGIT TYPES: " + str(diTypes))
print(len(FirstTimes))
evnum = evnums[ENTRYNUM]
binwidth = float(max(diT[ENTRYNUM]) - min(diT[ENTRYNUM]))/NBINS
#binedges = np.arange(min(diT[ENTRYNUM]), max(diT[ENTRYNUM]),binwidth) 
binedges = np.arange(7., 20.0, 13.0/NBINS)
#Make this into a numpy histogram and we do a fit to the histogram itself
timehist, binedges = np.histogram(FirstTimes,bins=binedges)
bincenters = []
for j,b in enumerate(binedges):
    if j == 0:
        continue
    bc = binedges[j] - (binedges[j] - binedges[j-1])
    bincenters.append(bc)
bincenters = np.array(bincenters)
print("HIST SUM: " + str(timehist.sum()))
popt, pcov = scp.curve_fit(thefunc, bincenters,timehist, p0=p0) 
print("BEST COEFF: " + str(popt))
#Now, get the best fit curve's points
bestfitfunc = thefunc(bincenters,popt[0], popt[1], popt[2],popt[3])
C2T = Chi2OverNDOF(bestfitfunc,timehist,timehist,NumVariables)
print("CHISQ TEST RESULT: " + str(C2T))
sns.set_style("whitegrid")
sns.axes_style("darkgrid")

xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
sns.set_palette(sns.xkcd_palette(xkcd_colors))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.hist(FirstTimes,bins=binedges, label="LAPPD Time Residuals")

plt.plot(bincenters,bestfitfunc,label=r'best fit, $\frac{\chi^{2}}{NDF}=%s$'%(str(np.round(C2T,2))))
#ax.plot(trueZvtx[0], trueXvtx[0], linestyle='none', 
#markersize=20, label="True Vertex", marker='*')
plt.ylabel("Entries")
plt.xlabel("Time Residual (ns)")
leg = ax.legend(loc=1)
leg.set_frame_on(True)
leg.draw_frame(True)
plt.title(r"Time Residual for event number %i"%(evnum)+"\n"+\
        "Landau distribution")
plt.show()
