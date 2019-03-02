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

def PMTHitCharge(Qs,typedata,title=None):
    binedges = np.arange(-0.0, 250.0, 50.0/NBINS)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    pmtind = np.where(typedata==0)[0]
    xkcd_colors = [ 'slate blue', 'warm pink', 'green', 'grass']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.hist(Qs[pmtind],bins=binedges, label="PMT Hit Charges")
    plt.ylabel("Entries")
    plt.xlabel("Hit charge")
    leg = ax.legend(loc=1)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    if title is None:
        plt.title("Hit charge distribution for PMT hits")
    else:
        plt.title(title)
    plt.show()

def LAPPDHitCharge(Qs,typedata,title=None):
    binedges = np.arange(0.0, 250.0, 50.0/NBINS)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    LAPPDind = np.where(typedata==1)[0]
    xkcd_colors = [ 'slate blue', 'warm pink', 'green', 'grass']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.hist(Qs[LAPPDind],bins=binedges, label="LAPPD Hit Charges")
    plt.ylabel("Entries")
    plt.xlabel("Hit charge")
    leg = ax.legend(loc=1)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    if title is None:
        plt.title("Hit charge distribution for LAPPD hits")
    else:
        plt.title(title)
    plt.show()

if __name__=='__main__':
    f = uproot.open("./ParaComb.root")
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
    PMTHitCharge(EventHitQs,EventdiType)
    LAPPDHitCharge(EventHitQs,EventdiType)
