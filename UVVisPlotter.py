#This script has functions needed to make a Cumulative Distribution Plot from
#Different variables output in  PhaseIITreeMaker root file.

import glob
import sys
import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_context('poster')
sns.set(font_scale=3.0)

def _GetSpec(filename):
    '''Given a filepath leading to a txt file with UV-Vis data, 
    collect all wavelengths and absorptions in an array.  return them.'''
    wvl=[]
    absorption=[]
    with open(filename,"r") as f:
        for line in f:
            if line.find("#") != -1:
                continue
            dataline = line.split(",")
            absnum = dataline[1].rstrip("\n")
            if absnum == "":
                print("No absorption value for this entry.  Ignoring")
                continue
            wvl.append(float(dataline[0]))
            absorption.append(float(absnum))
    return np.array(wvl),np.array(absorption)

if __name__=='__main__':
    print("USAGE: python UVVisPlotter.py [filepath] [plot_title]") 
    fpath = str(sys.argv[1])
    wvl,absorption = _GetSpec(fpath)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['slate blue', 'black', 'slate blue', 'warm pink', 'green', 'grass']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.axvline(300.0, color='black', linewidth=4, alpha=0.7)
    ax.plot(wvl,absorption, linewidth=4)
    plt.ylabel("Absorption")
    plt.xlabel("Wavelength (nm)")
    leg = ax.legend(loc=1,fontsize=14)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    plt.title(r" %s"%(str(sys.argv[2])))
    plt.show()
