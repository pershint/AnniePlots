import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.optimize as scp
import numpy as np
import sys
import copy

import lib.ROOTProcessor as rp
import SimpleKDE as skd

NBINS = 150

def GetLenQDistribution(f1data,indices=[]): 
    DiQ = np.array(f1data['digitQ'])
    DiType = np.array(f1data['digitType'])
    trueLen = np.array(f1data['trueTrackLengthInWater'])
    #Get the phi
    totalQs_pmt = [] 
    totalQs_lappd = []
    allQs_pmt = np.array([])
    allQs_lappd = np.array([])
    mu_length = [] 
    if len(indices)==0:
        indices = range(len(DiQ))
    for e in indices:
        thisQ = np.array(DiQ[e])
        typedata = np.array(DiType[e])
        #For each digit, we need the dir from trueVtx to the digit
        pmtind = np.where(typedata==0)[0]
        lappdind = np.where(typedata==1)[0]
        allQs_pmt = np.append(thisQ[pmtind],allQs_pmt)
        allQs_lappd = np.append(thisQ[lappdind],allQs_lappd)
        totalQs_pmt.append(np.sum(thisQ[pmtind]))
        totalQs_lappd.append(np.sum(thisQ[lappdind]))
        mu_length.append(trueLen[e])
    print("TOTAL CHARGE: " + str(np.array(totalQs_pmt)))
    print("TRACK LENGTHS: " + str(np.array(mu_length)))
    return allQs_pmt, allQs_lappd,np.array(totalQs_pmt),np.array(totalQs_lappd),np.array(mu_length)

def RunKDE(df,bandwidth=300,xb=100j,yb=100j):
    df_copy = copy.deepcopy(df)
    df_copy['PMT_sumQ'] = df_copy.PMT_sumQ/1000.
    MuEstimator = skd.KernelDensityEstimator(dataframe=df_copy)
    xx,yy,zz = MuEstimator.KDEEstimate2D(bandwidth,'muon_length','PMT_sumQ',xbins=xb,ybins=yb)
    zz=zz/np.max(zz)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors =  [ 'slate blue', 'green', 'grass','pink']
    sns.set_context("poster")
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    cbar = plt.contourf(xx,yy,zz,40,cmap='inferno')
    plt.colorbar()
    plt.title("Track length vs. total charge KDE for muon only, normalized to max of 1")
    plt.show()
    
    plt.show()
    return xx,yy,zz

if __name__ == '__main__':
    if str(sys.argv[1])=="--help":
        print(("USAGE: python HAP_new.py [variable_name_infile]"))
        sys.exit(0)
    f1 = str(sys.argv[1])

    #Process data into JSON objects
    mybranches = ['trueVtxX','trueVtxY','trueVtxZ','digitQ','digitType','trueTrackLengthInWater',
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
    mu_allQ_pmt, mu_allQ_lappd, mu_sumQ_pmt, mu_sumQ_lappd,mu_mulen = GetLenQDistribution(f1data,f1_mu_ind)
    pi_allQ_pmt, pi_allQ_lappd, pi_sumQ_pmt, pi_sumQ_lappd,pi_mulen = GetLenQDistribution(f1data,f1_pi_ind)
    #Apply a max charge cut; all events with Q>6000 considered muon + pion
    mu_sumQ_pmt_6000cutind = np.where(mu_sumQ_pmt<6000)[0]
    mu_sumQ_pmt = mu_sumQ_pmt[mu_sumQ_pmt_6000cutind]
    mu_mulen = mu_mulen[mu_sumQ_pmt_6000cutind]
    pi_sumQ_pmt_6000cutind = np.where(pi_sumQ_pmt<6000)[0]
    pi_sumQ_pmt = pi_sumQ_pmt[pi_sumQ_pmt_6000cutind]
    pi_mulen = pi_mulen[pi_sumQ_pmt_6000cutind]

    #Plot track length vs. sum of charge distributions 
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors =  [ 'grass', 'slate blue', 'green', 'pink']
    sns.set_context("poster")
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    data = {"PMT_sumQ":mu_sumQ_pmt, "muon_length":mu_mulen}
    mu_data = pd.DataFrame(data)
    g = sns.jointplot("muon_length","PMT_sumQ",data=mu_data,kind="hex",
            joint_kws=dict(gridsize=80),
            stat_func=None).set_axis_labels("Muon track length in tank (m)","PMT Total Charge (pe)")
    plt.subplots_adjust(left=0.2,right=0.8,
            top=0.85,bottom=0.2)
    cbar_ax = g.fig.add_axes([0.84,0.2,0.05,0.62])
    plt.colorbar(cax=cbar_ax)
    g.fig.suptitle("Total charge observed for a given muon track length \n" + \
            "Successful reconstruction, muon stops in MRD, muon only")
    plt.show()


    #Plot track length vs. sum of charge distributions 
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors =  [ 'grass', 'slate blue', 'green', 'pink']
    sns.set_context("poster")
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    data = {"PMT_sumQ":pi_sumQ_pmt, "muon_length":pi_mulen}
    pi_data = pd.DataFrame(data)
    g = sns.jointplot("muon_length","PMT_sumQ",data=pi_data,kind="hex",
            joint_kws=dict(gridsize=80),
            stat_func=None).set_axis_labels("Muon track length in tank (m)","PMT Total Charge (pe)")
    plt.subplots_adjust(left=0.2,right=0.8,
            top=0.85,bottom=0.2)
    cbar_ax = g.fig.add_axes([0.84,0.2,0.05,0.62])
    plt.colorbar(cax=cbar_ax)
    g.fig.suptitle("Total charge observed for a given muon track length \n" + \
            "Successful reconstruction, muon stops in MRD, muon + pions")
    plt.show()

#Let's use some KDEEEEEEEEEE
    xx,yy,zz = RunKDE(mu_data,bandwidth=10,xb=40j,yb=40j)
