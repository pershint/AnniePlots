#This script has functions needed to make a Cumulative Distribution Plot from
#Different variables output in  PhaseIITreeMaker root file.

import glob
import ROOT
import sys
import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_context('poster')

def _buildcsum(sorted_distribution):
    '''Given a distribution, returns the percentage from
    0 to 100 each data point is associated with'''
    c = np.arange(1, len(sorted_distribution) + 1, 1)
    h = c/(float(len(sorted_distribution)-1))
    return h*100.0

def _getCLdatavalue(CL,sorted_datadist):
    '''Given a fractional CL (ex: 0.683 for 68.3%), returns the
    value in the cumulative distribution's x-axis consistent with
    that CL'''
    CLindex = int(len(sorted_datadist)*CL)
    return sorted_datadist[CLindex]

if __name__=='__main__':
    print("USAGE: python CumulativeDistPlotter.py BRANCHNAME ROOTFILES") 
    f = uproot.open(str(sys.argv[2]))
    ftree = f.get("phaseII")
    ftree.items()
    branchname = sys.argv[1]
    branch = None
    try:
        branch = ftree.get(branchname)
        recoStatus = ftree.get("recoStatus")
    except:
        print("GIVE A VALID BRANCHNAME FOR THE THE INPUT ROOTFILE")
        print("USAGE: python CumulativeDistPlotter.py ROOTFILE BRANCHNAME")
        raise
    branch_arr_all = branch.array()
    status_arr = recoStatus.array()
    goodfits = np.where(status_arr == 0)[0]
    branch_arr = branch_arr_all[goodfits]
    sorted_brancharr = np.sort(branch_arr)
    csum = _buildcsum(sorted_brancharr)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['light eggplant', 'black', 'slate blue', 'warm pink', 'green', 'grass']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    CLvalue = _getCLdatavalue(0.68,sorted_brancharr)
    ax.plot(sorted_brancharr,csum, linewidth=6,label="MCTruth/Reco difference")
    ax.axvline(CLvalue, color='black', linewidth=6, label="68%% CL=%s"%(str(np.round(CLvalue,1))))
    plt.ylabel("Cumulative Distribution (%)")
    plt.xlabel(branchname + " (deg)")
    leg = ax.legend(loc=1,fontsize=14)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    plt.title("Cumulative distribution of variable %s,"%(branchname) + \
            " Full reconstruction chain \n 128 PMTs + 5 LAPPDs, recoStatus==0")
    plt.show()
    #Following would work if you want to combine individual ROOT files
    #in here, rather than using
    #$ hadd -k Newfile.root *.root
    #In the directory with your separate small ROOT files
    #print("USAGE: python CumulativeDistPlotter.py DIRECTORY")
    #filelist = glob.glob("%s*.root"%(direc))
    #direc = sys.argv[1:len(sys.argv)] 
#   #Let's try using the iterate method in uproot
#    for array in uproot.iterate("%s*.root"%(direc),"phaseII",["branchname1",
#        "branchname2"],entrysteps=10000):
#        print(array)
        #Need to figure out how to use the arrays though
