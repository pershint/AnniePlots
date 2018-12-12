import uproot
import HitChargePlotter as hcp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
sns.set_context('poster')
import sys
import numpy as np

NBINS=180


def PlotThetaDistribution(xdata, ydata, zdata,typedata,timedata,qdata):
    '''For the given digit data, plots the histogram of the 
    theta distribution for all hits'''
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
    y = np.array(y)
    theta = np.array(theta)
    pmtind = np.where(typedata==0)[0]
    LAPPDind = np.where(typedata==1)[0]
    time = np.array(timedata)
    charge = np.array(qdata)
    if len(LAPPDind) > 0:
        xkcd_colors = ['slate blue', 'warm pink', 'green', 'grass']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #plt.hist(LAPPD_SOL_wat_residuals,bins=lappd_binedges,label="n=1.33")
        thetabinedges = np.arange(-180.0, 180.0, 360.0/NBINS)
        plt.hist(theta[LAPPDind],bins=thetabinedges,label="LAPPD Thetas")
        plt.legend()
        plt.ylabel("Entries")
        plt.xlabel("LAPPD Thetas (degrees)")
        plt.title("Theta distribution of all LAPPD hits",fontsize=30)
        leg = ax.legend(loc=2)
        leg.set_frame_on(True)
        leg.draw_frame(True)
        plt.show()
    if len(pmtind) > 0:
        xkcd_colors = ['slate blue', 'warm pink', 'green', 'grass']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #plt.hist(LAPPD_SOL_wat_residuals,bins=lappd_binedges,label="n=1.33")
        thetabinedges = np.arange(-180.0, 180.0, 360.0/NBINS)
        plt.hist(theta,bins=thetabinedges,label="PMT Thetas")
        plt.legend()
        plt.ylabel("Entries")
        plt.xlabel("PMT Thetas (degrees)")
        plt.title("Theta distribution of all PMT hits with time residual < -1 ns",fontsize=30)
        leg = ax.legend(loc=2)
        leg.set_frame_on(True)
        leg.draw_frame(True)
        plt.show()


if __name__=='__main__':
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
    digitQ = ftree.get("digitQ")
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
    trueT = TrueVtxTime.array()
    print("ARRAY OF TRUE MUON VTX TIMES: " + str(trueT))
    diType = digitType.array()
    tot_entries = len(evn) 
    SOL_residuals = np.array([])
    types = np.array([]) 
    theX = np.array([])
    theY = np.array([]) 
    theZ = np.array([]) 
    theQ = np.array([])
    theT = np.array([])
    for e in xrange(tot_entries):
        #Want to calculate all hit residuals based on if
        #Light is moving at C or at C/1.33
        typedata = np.array(diType[e])
        pmtind = np.where(typedata==0)[0]
        lappdind = np.where(typedata==1)[0]
        this_TT = trueT[e]
        this_diQ = np.array(diQ[e])
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
        #pmt_res_SOL = res_SOL[pmtind]
        #PMT_SOL_residuals = np.append(pmt_res_SOL,PMT_SOL_residuals)
        #lappd_res_SOL = res_SOL[lappdind]
        #LAPPD_SOL_residuals = np.append(lappd_res_SOL,LAPPD_SOL_residuals)
        SOL_residuals = np.append(res_SOL, SOL_residuals)     
        types = np.append(typedata, types)
        theX = np.append(this_diX, theX)
        theY = np.append(this_diY, theY)
        theZ = np.append(this_diZ, theZ)
        theQ = np.append(this_diQ, theQ)
        theT = np.append(this_diT, theT)
    SOLB = np.where(SOL_residuals < -1.0)[0]
    PlotThetaDistribution(theX[SOLB], theY[SOLB], 
            theZ[SOLB],types[SOLB],theT[SOLB],
            theQ[SOLB])
    hcp.PMTHitCharge(theQ,types,title="Hit charge distribution for PMTs with time residual < -1 ns")
