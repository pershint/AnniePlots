import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
sns.set_context('poster')
import sys
import numpy as np

def IDsInRange(xdata, ydata, zdata, IDdata, ymin, ymax, thetamin, thetamax):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    theta = []
    for i in range(len(xdata)):
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
    if len(theta) == 0:
        print("No digits in range selected.  Continuing")
        return
    for j,ID in enumerate(IDdata):
        if ydata[j] >ymin and ydata[j] < ymax:
            if theta[j] > thetamin and theta[j] < thetamax:
                idCut.append(ID)
    plt.hist(idCut)
    ax.set_xlabel("LAPPD ID", fontsize=22)
    ax.set_ylabel("Digit counts", fontsize=22)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
    plt.title("LAPPDIDs in y and theta region identified")
    plt.show()

def XZ_ToTheta(xdata, zdata):
    theta = []
    for i in range(len(xdata)):
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

        theta.append(thistheta) #yeah that's confusing labeling
    return np.array(theta)

def YVSTheta(xdata, ydata, zdata,typedata,timedata,qdata):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    y = np.array(ydata)
    theta = XZ_ToTheta(xdata,zdata)
    pmtind = np.where(typedata==0)[0]
    LAPPDind = np.where(typedata==1)[0]
    time = np.array(timedata)
    charge = np.array(qdata)
    #scatter plot with time as color
    ourvmin = 2.0 #np.min(time)
    ourvmax = 11.0 #np.max(time)
    sc = None
    if len(LAPPDind) > 0:
        sc = plt.scatter(theta[LAPPDind],y[LAPPDind],c=time[LAPPDind], s=30,marker='o',label='LAPPD hit',cmap=cm.jet,vmin=ourvmin, vmax=ourvmax)
    if len(pmtind) > 0:
        sc = plt.scatter(theta[pmtind],y[pmtind],c=time[pmtind],s=200, marker='o',label='PMT hit',cmap=cm.jet,vmin=ourvmin, vmax=ourvmax)
    cbar = plt.colorbar(sc,label='Time (ns)')
    cbar.set_label(label='Time (ns)', size=30)
    cbar.ax.tick_params(labelsize=30)
    leg = ax.legend(loc=2)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    ax.set_xlabel("Theta (deg)", fontsize=34)
    ax.set_ylabel("Y (cm)", fontsize=34)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(30)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(30)
    plt.title("Hit times for a simulated muon production in ANNIE",fontsize=34)
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #scatter plot with charge as color
    ourvmin = np.min(charge)
    ourvmax = np.max(charge)
    if len(LAPPDind) > 0:
        sc = plt.scatter(theta[LAPPDind],y[LAPPDind],c=charge[LAPPDind], s=30, marker='o',label='LAPPD hit',cmap=cm.jet,vmin=ourvmin, vmax=ourvmax)
    if len(pmtind) > 0:
    
       sc = plt.scatter(theta[pmtind],y[pmtind],c=charge[pmtind],s=200, marker='o',label='PMT hit',cmap=cm.jet,vmin=ourvmin, vmax=ourvmax)
    cbar = plt.colorbar(sc,label='Charge')
    cbar.set_label(label='Charge', size=30)
    cbar.ax.tick_params(labelsize=30)
    leg = ax.legend(loc=2)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    ax.set_xlabel("Theta (deg)", fontsize=34)
    ax.set_ylabel("Y (cm)", fontsize=34)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
    plt.title("Hit charges for a simulated muon production in ANNIE",fontsize=34)
    plt.show()



if __name__=='__main__':
    f = uproot.open("./LongLAPPD_Default.root")
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
    digitTime = ftree.get("digitT")
    diTime = digitTime.array()
    digitQ = ftree.get("digitQ")
    diQ = digitQ.array()
    digitType = ftree.get("digitType")
    diType = digitType.array()
    digitID = ftree.get("digitDetID")
    diID = digitID.array()
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['eggplant', 'grass', 'vomit', 'warm pink', 'green', 'grass']
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

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(diX[eventnum],diY[eventnum], linestyle='none', marker='o',
    markersize=4, label="Digit Positions")
    ax.plot(trueXvtx[eventnum], trueYvtx[eventnum], linestyle='none', 
    markersize=20, label="True Vertex", marker='*')
    plt.ylabel("Y [cm]")
    plt.xlabel("X [cm]")
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
    
    YVSTheta(diX[eventnum],diY[eventnum],diZ[eventnum],diType[eventnum],diTime[eventnum],diQ[eventnum])
    IDsInRange(diX[eventnum],diY[eventnum],diZ[eventnum],diID[eventnum],-10,10,40, 60)
