#This script has functions needed to make a Cumulative Distribution Plot from
#Different variables output in  PhaseIITreeMaker root file.

import glob

import sys
import uproot
import lib.ROOTProcessor as rp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_context('poster')
sns.set(font_scale=4.0)

if __name__=='__main__':
    if str(sys.argv[1])=="--help":
        print(("USAGE: python LAPPDHitHists.py [filename] [title]"))
        sys.exit(0)
    f1 = str(sys.argv[1])
    title_string = str(sys.argv[2])
    #Process data into JSON objects
    mybranches = ['digitDetID','digitType','digitT']
    f1Processor = rp.ROOTProcessor(treename="phaseII")
    #branch = ['Pi0Count']
    f1Processor.addROOTFile(f1,branches_to_get=mybranches)
    f1data = f1Processor.getProcessedData()

    #Get our data arrays for making the nhit per LAPPD plot
    f1_IDs = np.array(f1data["digitDetID"])
    f1_Types = np.array(f1data["digitType"])
    f1_Times = np.array(f1data["digitT"])

    print(f1_Types[3])

    #Make a dictionary with the key as LAPPDID, value as nhits per event
    LAPPDHits = {}
    for event in range(len(f1_IDs)):
        #Get a dictionary with ID as key and nhits as value
        LAPPDHits_ThisEvent = {}
        LAPPDID_indices = np.where(np.array(f1_Types[event]) == 1)[0]
        LAPPDIDs = np.array(f1_IDs[event])
        IDs_unique = np.array(list(set(LAPPDIDs[LAPPDID_indices]))) #ugh
        for entry in f1_IDs[event]:
            if entry in IDs_unique and entry not in LAPPDHits_ThisEvent:
                LAPPDHits_ThisEvent[entry] = 1
            elif entry in IDs_unique and entry in LAPPDHits_ThisEvent:
                LAPPDHits_ThisEvent[entry]+=1
        #Now, we have to append the nhits to the full event dictionary
        print("LAPPDHITS_THISEVENT:")
        print(LAPPDHits_ThisEvent)
        for key in LAPPDHits_ThisEvent:
            if key in LAPPDHits:
                LAPPDHits[key].append(LAPPDHits_ThisEvent[key])
            else:
                LAPPDHits[key] = [LAPPDHits_ThisEvent[key]]

    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['purple', 'red',  'cobalt', 'blue', 'cobalt']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.xlabel('nhits')
    plt.ylabel('events')
    for key in LAPPDHits:
        if key not in [13,14,15]:
            continue
        #Plot a histogram of each LAPPD's nhits data
        plt.hist(LAPPDHits[key], range=(0,300), bins=150, 
                histtype='step',linewidth=4,label='LAPPD ID %s'%(str(key)))
    leg = ax.legend(loc=4,fontsize=24)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    plt.title(title_string)
    plt.show()
