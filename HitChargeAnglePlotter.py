import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.optimize as scp
import numpy as np
import sys

import lib.ROOTProcessor as rp
import SimpleKDE as skd

NBINS = 150

def GetPhiQDistribution(f1data,indices=[]): 
    DiZ = np.array(f1data['digitZ'])
    DiY = np.array(f1data['digitY'])
    DiX = np.array(f1data['digitX'])
    DiQ = np.array(f1data['digitQ'])
    DiType = np.array(f1data['digitType'])
    truDirZ = np.array(f1data['trueDirZ'])
    truDirY = np.array(f1data['trueDirY'])
    truDirX = np.array(f1data['trueDirX'])
    truX = np.array(f1data['trueVtxZ'])
    truY = np.array(f1data['trueVtxY'])
    truZ = np.array(f1data['trueVtxX'])
    #Get the phi
    allphis_pmt = np.array([])
    allQs_pmt = np.array([])
    allphis_lappd = np.array([])
    allQs_lappd = np.array([])
    if len(indices)==0:
        indices = range(len(DiZ))
    for e in indices:
        thisQ = np.array(DiQ[e])
        typedata = np.array(DiType[e])
        #For each digit, we need the dir from trueVtx to the digit
        thisDiX =np.array(DiX[e])
        thisDiY =np.array(DiY[e])
        thisDiZ =np.array(DiZ[e])
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
    return allphis_pmt, allQs_pmt, allphis_lappd, allQs_lappd

def RunKDE(df,bandwidth=300,xb=100j,yb=100j):
    MuEstimator = skd.KernelDensityEstimator(dataframe=df)
    xx,yy,zz = MuEstimator.KDEEstimate2D(bandwidth,'PMT_Phi','PMT_Charge',xbins=xb,ybins=yb)
    zz=zz/np.max(zz)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors =  [ 'slate blue', 'green', 'grass','pink']
    sns.set_context("poster")
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    cbar = plt.contourf(xx,yy,zz,40,cmap='inferno')
    plt.colorbar()
    plt.title("Hit charge vs. hit angle KDE for muon only, normalized to max of 1")
    plt.show()
    
    plt.show()
    return xx,yy,zz

if __name__ == '__main__':
    if str(sys.argv[1])=="--help":
        print(("USAGE: python HAP_new.py [variable_name_infile]"))
        sys.exit(0)
    f1 = str(sys.argv[1])

    #Process data into JSON objects
    mybranches = ['trueVtxX','trueVtxY','trueVtxZ','trueDirX','trueDirY','trueDirZ',
            'digitX','digitY','digitZ','digitQ','digitType',
            'recoVtxFOM','PiPlusCount','Pi0Count','PiMinusCount']
    f1Processor = rp.ROOTProcessor(treename="phaseII")
    #branch = ['Pi0Count']
    f1Processor.addROOTFile(f1,branches_to_get=mybranches)
    f1data = f1Processor.getProcessedData()

    #Split data into two subsets; with and without a pion
    f1_goodfitinFV = np.where((np.array(f1data['recoVtxFOM']) > 0) & 
            (np.array(f1data['trueVtxZ']) > 0) & (np.array(f1data['trueVtxY'])<50) & 
            (np.array(f1data['trueVtxY'])>-50))[0]
    f1_muononly = np.where((np.array(f1data['PiPlusCount']==0)) & 
            (np.array(f1data['Pi0Count'])==0) & 
            (np.array(f1data['PiMinusCount'])==0))[0]
    f1_pions = np.where((np.array(f1data['PiPlusCount'])>0) | 
            (np.array(f1data['Pi0Count'])>0) | 
            (np.array(f1data['PiMinusCount'])>0))[0]
    f1_mu_ind = np.intersect1d(f1_goodfitinFV,f1_muononly) 
    f1_pi_ind = np.intersect1d(f1_goodfitinFV,f1_pions) 

    #Get the charge vs. phi data for both total datasets
    mu_phi_pmt, mu_Q_pmt, mu_phi_lappd, mu_Q_lappd = GetPhiQDistribution(f1data,f1_mu_ind)
    pi_phi_pmt, pi_Q_pmt, pi_phi_lappd, pi_Q_lappd = GetPhiQDistribution(f1data,f1_pi_ind)

    #Apply a max charge cut; all PMTs with Q>500 considered muon + pion light
    mu_Q_pmt_500cutind = np.where(mu_Q_pmt<500)[0]
    mu_phi_pmt = mu_phi_pmt[mu_Q_pmt_500cutind]
    mu_Q_pmt = mu_Q_pmt[mu_Q_pmt_500cutind]
    pi_Q_pmt_500cutind = np.where(pi_Q_pmt<500)[0]
    pi_phi_pmt = pi_phi_pmt[pi_Q_pmt_500cutind]
    pi_Q_pmt = pi_Q_pmt[pi_Q_pmt_500cutind]


    #Plot Q vs. Phi distributions 
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors =  ['slate blue', 'green', 'grass','pink']
    sns.set_context("poster")
    sns.set_palette(sns.xkcd_palette(xkcd_colors))

    data = {"PMT_Phi":mu_phi_pmt, "PMT_Charge":mu_Q_pmt}
    mu_data = pd.DataFrame(data)
    g = sns.jointplot("PMT_Phi","PMT_Charge",data=mu_data,kind="hex",
            joint_kws=dict(gridsize=80),
            stat_func=None).set_axis_labels("PMT Phi from vertex dir. (deg)","PMT Charge (pe)")
    plt.subplots_adjust(left=0.2,right=0.8,
            top=0.85,bottom=0.2)
    cbar_ax = g.fig.add_axes([0.84,0.2,0.05,0.62])
    plt.colorbar(cax=cbar_ax)
    g.fig.suptitle("Distribution of PMT charges relative to muon direction \n" + \
            "Successful reconstruction, muon stops in MRD, muon only")
    plt.show()

#Plot Q vs. Phi distributions 
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors =  [ 'slate blue', 'green', 'grass','pink']
    sns.set_context("poster")
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    data = {"PMT_Phi":pi_phi_pmt, "PMT_Charge":pi_Q_pmt}
    pi_data = pd.DataFrame(data)
    g = sns.jointplot("PMT_Phi","PMT_Charge",data=pi_data,kind="hex",
            joint_kws=dict(gridsize=80),
            stat_func=None).set_axis_labels("PMT Phi from vertex dir. (deg)","PMT Charge (pe)")
    plt.subplots_adjust(left=0.2,right=0.8,
            top=0.85,bottom=0.2)
    cbar_ax = g.fig.add_axes([0.84,0.2,0.05,0.62])
    plt.colorbar(cax=cbar_ax)
    g.fig.suptitle("Distribution of PMT charges relative to muon direction \n" + \
            "Successful reconstruction, muon stops in MRD, muon + pions")
    plt.show()

#Let's use some KDEEEEEEEEEE
    xx,yy,zz = RunKDE(mu_data,bandwidth=10,xb=40j,yb=40j)
