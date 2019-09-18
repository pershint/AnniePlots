#This script has functions needed to make a Cumulative Distribution Plot from
#Different variables output in  PhaseIITreeMaker root file.

import glob

import sys
import uproot
import lib.ROOTProcessor as rp
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import numpy as np
import pandas as pd

sns.set_context('poster')
sns.set(font_scale=2.0)

PMTRADIUS = 100.
PMTHALFHEIGHT = 145.

def CartToPolar(X,Z):
    '''
    Given the data, produce radial and angular bin values where the
    value for each bin is the mean of the given variable in the valid range.
    '''
    r = []
    phi = []
    for i in range(len(X)):
        thisphi = None
        x = X[i]
        z = Z[i]
        isnegative = False
        if x <0:
            isnegative = True
            x = -x
        rad = np.sqrt(z**2 + x**2)
        r.append(rad)
        if rad == 0:
            thisphi = np.arccos(0)*180.0/np.pi
        else:
            thisphi = np.arccos(z/rad)*180.0/np.pi 
        if isnegative:
            thisphi = (360.0 - thisphi)
            #thisphi = (180.0 + thisphi)
        phi.append(thisphi)
    return np.array(r),np.array(phi)

def XZ_means_radial(X,Z,variable,ang_bins=10,rad_bins=5):
    '''
    Given the data, produce radial and angular bin values where the
    value for each bin is the mean of the given variable in the valid range.
    '''
    r_avgbin = []
    phi_avgbin = []
    var_avgbin = []
    r,phi = CartToPolar(X,Z)
    angular_bins = np.arange(0,360.,360./ang_bins)
    print("ANGULAR BINS ARE: " + str(angular_bins))
    ang_res = 360./ang_bins
    radial_bins = np.arange(0,PMTRADIUS,PMTRADIUS/rad_bins)
    rad_res = PMTRADIUS/rad_bins
    for i,aval in enumerate(angular_bins):
        if i==0: continue
        for j,rval in enumerate(radial_bins):
            if j==0: continue
            r_avgbin.append(radial_bins[j-1] + rad_res/2.)
            phi_avgbin.append(angular_bins[i-1] + ang_res/2.)
            #Now, Get the variables fitting into this r/angle bin
            binvalinds1 = np.where((r > radial_bins[j-1]) & (r < radial_bins[j]))[0]
            binvalinds2 = np.where((phi > angular_bins[i-1]) & (phi < angular_bins[i]))[0]
            thebinvals = np.intersect1d(binvalinds1,binvalinds2)
            theavgval =np.average(variable[thebinvals])
            var_avgbin.append(theavgval)
    return np.array(r_avgbin),np.array(phi_avgbin),np.array(var_avgbin)

def XZ_means(X,Z,variable,x_bins=10,z_bins=5):
    '''
    Given the data, produce radial and angular bin values where the
    value for each bin is the mean of the given variable in the valid range.
    '''
    x_bins = []
    z_bins = []
    for i,aval in enumerate(angular_bins):
        if i==0: continue
        for j,rval in enumerate(radial_bins):
            if j==0: continue
            r_avgbin.append(radial_bins[j-1] + rad_res/2.)
            phi_avgbin.append(angular_bins[i-1] + ang_res/2.)
            #Now, Get the variables fitting into this r/angle bin
            binvalinds1 = np.where((r > radial_bins[j-1]) & (r < radial_bins[j]))[0]
            binvalinds2 = np.where((phi > angular_bins[i-1]) & (phi < angular_bins[i]))[0]
            thebinvals = np.intersect1d(binvalinds1,binvalinds2)
            theavgval =np.average(variable[thebinvals])
            var_avgbin.append(theavgval)
    return np.array(r_avgbin),np.array(phi_avgbin),np.array(var_avgbin)


def ContourMap_XZSlice(X,Y,Z,deltaVtxR,yrange=[-50.0,50.0],ngridx=60, ngridz=60):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    x, z, dvr = [], [], []
    valid_ind = np.where((Y>yrange[0]) & (Y<yrange[1]))[0]
    x = np.array(X[valid_ind])
    z = np.array(Z[valid_ind])
    dvr = np.array(deltaVtxR[valid_ind])
    #Perform linear interpolation of the data (x,z) on a grid
    #Defined by (xi, zi) as seen in the Matplotlib examples
    triang = tri.Triangulation(z,x)
    interpolator = tri.LinearTriInterpolator(triang, dvr)
    xi = np.linspace(2*np.min(x), 2*np.max(x), ngridx)
    zi = np.linspace(2*np.min(z), 2*np.max(z), ngridz)
    Zi, Xi = np.meshgrid(zi,xi)
    dvri = interpolator(Zi, Xi)
    print(dvri)
    dvri = gaussian_filter(dvri, sigma=0.8)
    print(dvri)
    contour1 = ax1.contourf(zi, xi, dvri, 10, cmap=plt.cm.jet)
    fig.colorbar(contour1, ax=ax1)
    ax1.set_title("Best fit contour using linear interpolation")
    plt.show()

def GridMap_XZSlice(X,Y,Z,deltaVtxR,yrange=[-50.0,50.0],ngridx=30, ngridz=30):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    x, z, dvr = [], [], []
    valid_ind = np.where((Y>yrange[0]) & (Y<yrange[1]))[0]
    x = np.array(X[valid_ind])
    z = np.array(Z[valid_ind])
    dvr = np.array(deltaVtxR[valid_ind])
    #Perform linear interpolation of the data (x,z) on a grid
    #Defined by (xi, zi) as seen in the Matplotlib examples
    xi = np.linspace(np.min(x), np.max(x), ngridx)
    zi = np.linspace(np.min(z), np.max(z), ngridz)
    print("XI: " + str(xi))
    print("ZI: " + str(zi))
    Zi, Xi = np.meshgrid(zi,xi)
    points = (z,x)
    grid_test = griddata(points, dvr, (Zi,Xi), method='linear')
    grid_test = gaussian_filter(grid_test,sigma=0.8)
    im = plt.imshow(grid_test, extent = (np.min(z), np.max(z), np.min(x), np.max(x)))
    cbar = fig.colorbar(im)
    cbar.set_label("Reco. resolution (cm)")
    ax1.set_title("Reco. resolution in tank, yrange [-50,50] cm\n(fit to data using linear interpolation)")
    ax1.set_ylabel("True X-Position (cm)")
    ax1.set_xlabel("True Z-Position (cm)")
    plt.show()

def FOMMap_XZSlice(X,Y,Z,FOM,yrange=[-50.0,50.0],ngridx=30, ngridz=30):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    x, z, dvr = [], [], []
    valid_ind = np.where((Y>yrange[0]) & (Y<yrange[1]))[0]
    x = np.array(X[valid_ind])
    z = np.array(Z[valid_ind])
    fom = np.array(FOM[valid_ind])
    #Perform linear interpolation of the data (x,z) on a grid
    #Defined by (xi, zi) as seen in the Matplotlib examples
    xi = np.linspace(np.min(x), np.max(x), ngridx)
    zi = np.linspace(np.min(z), np.max(z), ngridz)
    Zi, Xi = np.meshgrid(zi,xi)
    points = (z,x)
    grid_test = griddata(points, fom, (Zi,Xi), method='linear')
    grid_test = gaussian_filter(grid_test,sigma=0.8)
    im = plt.imshow(grid_test, extent = (np.min(z), np.max(z), np.min(x), np.max(x)))
    cbar = fig.colorbar(im)
    cbar.set_label("Reco. resolution (cm)")
    ax1.set_title("Reco. FOM in tank, yrange [-50,50] cm\n(fit to data using linear interpolation)")
    ax1.set_ylabel("True X-Position (cm)")
    ax1.set_xlabel("True Z-Position (cm)")
    plt.show()

def TwoDHist_XZSlice(X,Y,Z,yrange=[-30.0,30.0]):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    x, z = [], []
    valid_ind = np.where((Y>yrange[0]) & (Y<yrange[1]))[0]
    x = np.array(X[valid_ind])
    z = np.array(Z[valid_ind])
    data = {"true_X":x, "true_Z":z}
    pd_data = pd.DataFrame(data)
    g = sns.jointplot("true_Z","true_X",data=pd_data,kind="hex",
            joint_kws=dict(gridsize=15),
            stat_func=None).set_axis_labels("True Z-Position (cm)","True X-Position (cm)")
    plt.subplots_adjust(left=0.2,right=0.8,
            top=0.95,bottom=0.2)
    cbar_ax = g.fig.add_axes([0.85,0.2,0.05,0.62])
    plt.colorbar(cax=cbar_ax)
    g.fig.suptitle("Histogram of events in Y-range[%d,%d]\n"%(yrange[0],yrange[1]) + \
            "Successful reconstruction, muon stops in MRD")
    plt.show()

if __name__=='__main__':
    if str(sys.argv[1])=="--help":
        print("USAGE: python RecoEfficiencyPlots.py [file1.root] ")
        sys.exit(0)
    f1 = str(sys.argv[1])
    #Process data into JSON objects
    mybranches = ['recoVtxFOM','recoVtxX','recoVtxY','recoVtxZ',
            'trueVtxX','trueVtxY','trueVtxZ','deltaAngle',
            'deltaVtxR']
    f1Processor = rp.ROOTProcessor(treename="phaseII")
    #branch = ['Pi0Count']
    f1Processor.addROOTFile(f1,branches_to_get=mybranches)
    f1data = f1Processor.getProcessedData()
    f1_fom = np.array(f1data['recoVtxFOM'])
    #Get only good fits
    goodfit_inds = np.where(f1_fom>0)[0]
    f1_recoVtxX = np.array(f1data['recoVtxX'])[goodfit_inds]
    f1_recoVtxY = np.array(f1data['recoVtxY'])[goodfit_inds]
    f1_recoVtxZ = np.array(f1data['recoVtxZ'])[goodfit_inds]
    f1_trueVtxX = np.array(f1data['trueVtxX'])[goodfit_inds]
    f1_trueVtxY = np.array(f1data['trueVtxY'])[goodfit_inds]
    f1_trueVtxZ = np.array(f1data['trueVtxZ'])[goodfit_inds]
    f1_fom = f1_fom[goodfit_inds]
    f1_deltaangle = np.array(f1data['deltaAngle'])[goodfit_inds]
    f1_deltar = np.array(f1data['deltaVtxR'])[goodfit_inds]
    print("NUM EVENTS AFTER ONLY GETTING GOOD FITS: %i"%(len(f1_deltar))) 
    
    rbin, phibin, valavg = XZ_means_radial(f1_trueVtxX,f1_trueVtxZ,f1_deltar,ang_bins=30,rad_bins=10)
    zbin = np.array(rbin)*np.cos(np.array(phibin))
    xbin = np.array(rbin)*np.sin(np.array(phibin))
    ybin = np.zeros(len(zbin))
    ContourMap_XZSlice(f1_trueVtxX,f1_trueVtxY,f1_trueVtxZ,f1_deltar,yrange=[-50,50])
    GridMap_XZSlice(f1_trueVtxX,f1_trueVtxY,f1_trueVtxZ,f1_deltar,yrange=[-50,50])
    FOMMap_XZSlice(f1_trueVtxX,f1_trueVtxY,f1_trueVtxZ,f1_fom,yrange=[-50,50])
    TwoDHist_XZSlice(f1_trueVtxX,f1_trueVtxY,f1_trueVtxZ,yrange=[-50,50])
