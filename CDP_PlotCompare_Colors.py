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
sns.set(font_scale=3.0)

def _rebin(cumulative_distribution,entrys_per_value):
    '''Given an array, take the mean of #entrys_per_value given and
    set that as the new bin value'''
    csum = np.array(cumulative_distribution)
    csum_rebin = []
    for j in np.arange(0,len(cumulative_distribution),entrys_per_value):
        if j == 0:
            continue
        newbin = np.average(csum[j-1-entrys_per_value:j-1])
        csum_rebin.append(newbin)
    return np.array(csum_rebin)


def _buildcsum(sorted_distribution):
    '''Given a distribution, returns the percentage from
    0 to 100 each data point is associated with'''
    c = np.arange(1, len(sorted_distribution) + 1, 1)
    h = c/(float(len(sorted_distribution)-1))
    return h*100.0

def _buildcsumbins(sorted_distribution,dist_per_bin,binrange=None):
    '''Given a distribution, returns the percentage from
    0 to 100 each data point is associated with'''
    bin_height = []
    if binrange is None:
        bin_lefts = np.arange(np.min(sorted_distribution),np.max(sorted_distribution),dist_per_bin)
    else:
        bin_lefts = np.arange(binrange[0], binrange[1]+dist_per_bin,dist_per_bin)
    bin_rights = bin_lefts + dist_per_bin
    for j,val in enumerate(bin_lefts):
        bin_right = val + dist_per_bin
        in_bin_inds = np.where((sorted_distribution > val) & (sorted_distribution<bin_right))[0]
        bin_frac = float(len(in_bin_inds))/float(len(sorted_distribution))
        if j==0:
            bin_height.append(bin_frac)
        else:
            bin_height.append(bin_frac + bin_height[j-1])
    bin_height = np.array(bin_height)*100.0
    return bin_height, bin_lefts, bin_rights

def _buildbins(sorted_distribution,dist_per_bin,binrange=None):
    '''Given a distribution, returns the percentage from
    0 to 100 each data point is associated with'''
    bin_height = []
    if binrange is None:
        bin_lefts = np.arange(np.min(sorted_distribution),np.max(sorted_distribution),dist_per_bin)
    else:
        bin_lefts = np.arange(binrange[0], binrange[1]+dist_per_bin,dist_per_bin)
    bin_rights = bin_lefts + dist_per_bin
    for j,val in enumerate(bin_lefts):
        bin_right = val + dist_per_bin
        in_bin_inds = np.where((sorted_distribution > val) & (sorted_distribution<bin_right))[0]
        bin_val = float(len(in_bin_inds))
        bin_height.append(bin_val)
    bin_height = np.array(bin_height)
    return bin_height, bin_lefts, bin_rights

def _getCLdatavalue(CL,sorted_datadist):
    '''Given a fractional CL (ex: 0.683 for 68.3%), returns the
    value in the cumulative distribution's x-axis consistent with
    that CL'''
    CLindex = int(len(sorted_datadist)*CL)
    return sorted_datadist[CLindex]

def _getCLdatavaluebin(CL,CLbins,bin_lefts):
    '''Given a fractional CL (ex: 0.683 for 68.3%), returns the
    value in the cumulative distribution's x-axis consistent with
    that CL'''
    thewin=np.where(CLbins>CL)[0][0]
    return bin_lefts[thewin], CLbins[thewin]

if __name__=='__main__':
    if str(sys.argv[1])=="--help":
        print(("USAGE: python CDP_PlotCompare.py [variable_name_infile] [variable_unit] [file1.root] "+
           "[file1_legend] [file2.root] [file2_legend]"))
        sys.exit(0)
    f1 = str(sys.argv[3])
    f1_name = str(sys.argv[4])
    f2 = str(sys.argv[5])
    f2_name = str(sys.argv[6])
    var_unit = str(sys.argv[2])
    #Process data into JSON objects
    mybranches = ['%s'%(str(sys.argv[1])),'recoVtxFOM','recoDirZ']
    f1Processor = rp.ROOTProcessor(treename="phaseII")
    #branch = ['Pi0Count']
    f1Processor.addROOTFile(f1,branches_to_get=mybranches)
    f1data = f1Processor.getProcessedData()
    print(f1data)

    f2Processor = rp.ROOTProcessor(treename="phaseII")
    #branch = ['Pi0Count']
    f2Processor.addROOTFile(f2,branches_to_get=mybranches)
    f2data = f2Processor.getProcessedData()

    f1_data = np.array(f1data[str(sys.argv[1])])
    f2_data = np.array(f2data[str(sys.argv[1])])
    f1_dirZ_arr = np.array(f1data['recoDirZ'])
    f2_dirZ_arr = np.array(f2data['recoDirZ'])
    f1_fom_arr = np.array(f1data['recoVtxFOM'])
    f2_fom_arr = np.array(f2data['recoVtxFOM'])
    f1_goodfitind = np.where((f1_fom_arr > 0) & (f1_dirZ_arr > 0))[0]
    f2_goodfitind = np.where((f2_fom_arr > 0) & (f2_dirZ_arr > 0))[0]
    f1_goodfits = f1_data[f1_goodfitind]
    f2_goodfits = f2_data[f2_goodfitind]
    sorted_f1 = np.sort(f1_goodfits)
    sorted_f2 = np.sort(f2_goodfits)
    binwidth=0
    themax =0
    if str(sys.argv[1]) == 'deltaVtxR':
        binwidth=1
    if str(sys.argv[1]) == 'deltaAngle':
        binwidth=1
    if str(sys.argv[1]) == 'deltaVtxR':
        variable_str = r'$\Delta$r'
        themax = 100
    if str(sys.argv[1]) == 'deltaAngle':
        variable_str = r'$\Delta \phi$'
        themax = 50
    f1_csum,f1_binleft,f1_binright = _buildcsumbins(sorted_f1,binwidth,binrange=[0,themax])
    print("f1_csum: " + str(f1_csum))
    print("f1_binleft: " + str(f1_binleft))
    print("f1_csum: " + str(f1_csum))
    print("f1_binleft: " + str(f1_binleft))
    f2_csum,f2_binleft,f2_binright = _buildcsumbins(sorted_f2,binwidth,binrange=[0,themax])
    f1_CLvalue, f1_CLvalueheight = _getCLdatavaluebin(68,f1_csum,f1_binleft)
    f2_CLvalue, f2_CLvalueheight = _getCLdatavaluebin(68,f2_csum,f2_binleft)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['dark teal' for x in range(len(f1_csum)*2)]
    csum2_colors = ['adobe' for x in range(len(f2_csum)*2)]
    colors = xkcd_colors + csum2_colors
    sns.set_palette(sns.xkcd_palette(colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    maxclwidth = np.max([f1_CLvalue, f2_CLvalue])
    ax.axhline(68, xmin=0, xmax=(maxclwidth/themax),color='black', linestyle='--',linewidth=4, label="68% CL",alpha=0.75)
    ax.axvline(f1_CLvalue, ymin=0, ymax=f1_CLvalueheight/100.,color='black', linestyle='--',linewidth=4, alpha=0.75)
    ax.axvline(f2_CLvalue, ymin=0, ymax=f2_CLvalueheight/100.,color='black', linestyle='--',linewidth=4, alpha=0.75)
    for j,val in enumerate(f1_csum):
        if j == len(f1_csum)-1:
            break
        elif j == 0:
            ax.plot([0,f1_binleft[j]],[0,0],linewidth=6,linestyle='-',label=f1_name)
            ax.plot([f1_binleft[j],f1_binleft[j]],[0,val],linewidth=6,linestyle='-')
            ax.plot([f1_binleft[j],f1_binright[j]],[val,val],linewidth=6,linestyle='-',label=f1_name)
            ax.plot([f1_binright[j],f1_binright[j]],[val,f1_csum[j+1]],linewidth=6,linestyle='-')
        else:
            ax.plot([f1_binleft[j],f1_binright[j]],[val,val],linewidth=6,linestyle='-')
            ax.plot([f1_binright[j],f1_binright[j]],[val,f1_csum[j+1]],linewidth=6,linestyle='-')
    for j,val in enumerate(f2_csum):
        if j == len(f2_csum)-1:
            break
        elif j == 0:
            ax.plot([0,f2_binleft[j]],[0,0],linewidth=6,linestyle='-',label=f2_name)
            ax.plot([f2_binleft[j],f2_binleft[j]],[0,val],linewidth=6,linestyle='-')
            ax.plot([f2_binleft[j],f2_binright[j]],[val,val],linewidth=6,linestyle='-',label=f2_name)
            ax.plot([f2_binright[j],f2_binright[j]],[val,f2_csum[j+1]],linewidth=6,linestyle='-')
        else:
            ax.plot([f2_binleft[j],f2_binright[j]],[val,val],linewidth=6,linestyle='-')
            ax.plot([f2_binright[j],f2_binright[j]],[val,f2_csum[j+1]],linewidth=6,linestyle='-')
    plt.ylabel("%% Probability/%s %s "%(str(binwidth),str(sys.argv[2])))
    plt.xlim(0,themax)
    plt.xlabel(variable_str+" [%s]"%(var_unit))
    plt.ylim(0,100)
    leg = ax.legend(loc=4,fontsize=24)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    plt.title("Reconstruction resolution of variable %s"%(variable_str))
    plt.show()
    
    #Make data histogram plot
    if str(sys.argv[1]) == 'deltaVtxR':
        variable_str = r'$\Delta$r'
        themax = 100
    if str(sys.argv[1]) == 'deltaAngle':
        variable_str = r'$\Delta \phi$'
        themax = 50
    f2_hist,f2_binleft,f2_binright = _buildbins(sorted_f2,binwidth,binrange=[0,themax])
    f1_hist,f1_binleft,f1_binright = _buildbins(sorted_f1,binwidth,binrange=[0,themax])
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['dark teal' for x in range(len(f1_hist)*2)]
    csum2_colors = ['adobe' for x in range(len(f2_hist)*2)]
    colors = xkcd_colors + csum2_colors
    sns.set_palette(sns.xkcd_palette(colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for j,val in enumerate(f1_hist):
        if j == len(f1_hist)-1:
            break
        elif j == 0:
            ax.plot([0,f1_binleft[j]],[0,0],linewidth=6,linestyle='-',label=f1_name)
            ax.plot([f1_binleft[j],f1_binleft[j]],[0,val],linewidth=6,linestyle='-')
            ax.plot([f1_binleft[j],f1_binright[j]],[val,val],linewidth=6,linestyle='-',label=f1_name)
            ax.plot([f1_binright[j],f1_binright[j]],[val,f1_hist[j+1]],linewidth=6,linestyle='-')
        else:
            ax.plot([f1_binleft[j],f1_binright[j]],[val,val],linewidth=6,linestyle='-')
            ax.plot([f1_binright[j],f1_binright[j]],[val,f1_hist[j+1]],linewidth=6,linestyle='-')
    for j,val in enumerate(f2_hist):
        if j == len(f2_hist)-1:
            break
        elif j == 0:
            ax.plot([0,f2_binleft[j]],[0,0],linewidth=6,linestyle='-',label=f2_name)
            ax.plot([f2_binleft[j],f2_binleft[j]],[0,val],linewidth=6,linestyle='-')
            ax.plot([f2_binleft[j],f2_binright[j]],[val,val],linewidth=6,linestyle='-',label=f2_name)
            ax.plot([f2_binright[j],f2_binright[j]],[val,f2_hist[j+1]],linewidth=6,linestyle='-')
        else:
            ax.plot([f2_binleft[j],f2_binright[j]],[val,val],linewidth=6,linestyle='-')
            ax.plot([f2_binright[j],f2_binright[j]],[val,f2_hist[j+1]],linewidth=6,linestyle='-')
    plt.ylabel("Counts/%s %s "%(str(binwidth),str(sys.argv[2])))
    plt.xlim(0,themax)
    plt.xlabel(variable_str+" [%s]"%(var_unit))
    themax = np.max([np.max(f1_hist),np.max(f2_hist)])
    plt.ylim(0,themax+10)
    leg = ax.legend(loc=1,fontsize=24)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    plt.title("Counts per bin for reconstructed variable %s"%(variable_str))
    plt.show()
