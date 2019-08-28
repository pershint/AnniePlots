import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import pandas as pd

import lib.ROOTProcessor as rp

sns.set_context('poster')
import sys
import numpy as np

def XZ_ToTheta(xdata, zdata):
    theta = []
    for i in range(len(xdata)):
        thistheta = None
        x = xdata[i]
        z = zdata[i]
        isnegative = False
        print("THIS X,Z: %d,%d"%(x,z))
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

        theta.append(thistheta) #yeah that's confusing labeling
        print("THIS THETA: " + str(thistheta))
    return np.array(theta)

def YVSTheta_all(xdata, ydata, zdata,typedata,iddata):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ids_plotted = []
    all_xs = []
    all_zs = []
    all_ys = []
    all_thetas = []
    all_rs = []
    print(len(xdata))
    for ev in range(len(xdata)):
        evx = xdata[ev]
        evy = ydata[ev]
        evz = zdata[ev]
        evtype = typedata[ev]
        evid = iddata[ev]
        for j,hit in enumerate(evid):
            if hit in ids_plotted or evtype[j]==1: 
                continue
            if abs(evy[j]) > 130:
                continue
            else:
                ids_plotted.append(hit)
                theta = XZ_ToTheta([evx[j]],[evz[j]])
                thisr = np.sqrt(evx[j]**2 + evz[j]**2)
                all_ys.append(evy[j])
                all_zs.append(evz[j])
                all_xs.append(evx[j])
                all_thetas.append(theta[0])
    print("ALL TANK WALL IDS: ")
    print(ids_plotted)
    #Let's cheese in that one missing PMT... weird...
    ids_plotted.append(-1)
    all_ys.append(103)
    all_thetas.append(-167)
    sc = plt.scatter(all_thetas,all_ys,c=ids_plotted, s=200, marker='o',label='PMT positions',
            vmin=0, vmax= 100)
    cbar = plt.colorbar(sc,label='PMT ID')
    cbar.set_label(label='PMT ID', size=30)
    cbar.ax.tick_params(labelsize=30)
    leg = ax.legend(loc=2)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    ax.set_xlabel("Theta (deg)", fontsize=34)
    ax.set_ylabel("Y (cm)", fontsize=34)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(30)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(30)
    plt.title("Positions of all PMTs hit somewhere in this file",fontsize=34)
    plt.show()
    plt.scatter(all_zs,all_xs,s=200, marker='o',label='PMT positions')
    leg = ax.legend(loc=2)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    ax.set_xlabel("Z (cm)", fontsize=34)
    ax.set_ylabel("X (cm)", fontsize=34)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(30)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(30)
    plt.title("Positions of all PMTs hit somewhere in this file",fontsize=34)
    plt.show()
    return np.array(ids_plotted), np.array(all_ys), np.array(all_thetas)

def PositionsToPixels(ids, ys, thetas):
    Y_ranges = [[90,120], [50,80], [10,40], [-40,-10], [-80, -50], [-120,-90]]

    Th_ranges = [ [-180,-160],[-160,-135], [-135,-110],  [-110, -85],
            [-85, -65], [-65, -50], [-50, -20], [-20,0], [0,20], 
             [20, 50],[50, 65], [65, 85], [85, 110], [110, 135], 
             [135, 160],[160, 180]]
    #Now, we need to map these IDs to pixel values
    #Start by sorting everything in terms of smallest ys
    print("SORTED THETAS")
    print(sorted(thetas))
    pixel_map = {"xpixel": [], "ypixel": [], "id": []}
    for j,yrange in enumerate(Y_ranges):
        for k,thrange in enumerate(Th_ranges):
            print("THE YRANGE: %s"%(str(yrange)))
            print("THE THRANGE: %s"%(str(thrange)))
            ranged_yind = np.where((ys>yrange[0]) & (ys<yrange[1]))[0]
            ranged_thind = np.where((thetas>thrange[0]) & (thetas<thrange[1]))[0]
            print("INDS IN YRANGE: %s"%(str(ranged_yind)))
            print("INDS IN THRANGE: %s"%(str(ranged_thind)))
            theindex = np.intersect1d(ranged_yind,ranged_thind)
            if len(theindex) > 1:
                print("OH CRAP, YOUR Y-THETA RANGE HAS MORE THAN ONE ENTRY...")
            pixel_map["id"].append(ids[theindex[0]])
            pixel_map["xpixel"].append(k)
            pixel_map["ypixel"].append(j)
    pixel_map_pandas = pd.DataFrame(pixel_map)
    print(pixel_map_pandas)
    pmp = pixel_map_pandas.pivot(index="ypixel",columns="xpixel",values="id")
    ax = sns.heatmap(pmp)
    plt.show()
    return pixel_map




if __name__=='__main__':

    if str(sys.argv[1])=="--help":
        print("USAGE: python PMTMap.py [file1.root] ")
        sys.exit(0)
    f1 = str(sys.argv[1])
    #Process data into JSON objects
    mybranches = ['digitX','digitY','digitZ',
            'digitType','digitDetID']
    f1Processor = rp.ROOTProcessor(treename="phaseII")
    #branch = ['Pi0Count']
    f1Processor.addROOTFile(f1,branches_to_get=mybranches)
    f1data = f1Processor.getProcessedData()
    diX = np.array(f1data['digitX'])
    diY = np.array(f1data['digitY'])
    diZ = np.array(f1data['digitZ'])
    diType = np.array(f1data['digitType'])
    diID = np.array(f1data['digitDetID'])

    ids, ys, thetas = YVSTheta_all(diX,diY,diZ,diType,diID)
    #Now, create the map of IDs to pixel indices
    pixel_map = PositionsToPixels(ids, ys, thetas)
    print(pixel_map)
    #with open("./pixel_map.json","w") as f:
    #    json.dump(f,pixel_map,indent=4,sort_keys=True)
