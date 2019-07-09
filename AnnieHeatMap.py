#When run in a directory containing all timebin results (outputs from 
#main.py), produces a table showing the data cleaning sacrifice and
#Data/MC ratio for each timebin.  

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import HitsPlotter as hp
sns.set_context('poster')
sns.set_style("darkgrid")
import os,sys,time
import glob
import json
import pandas as pd
import uproot

basepath = os.path.dirname(__file__)
MAINDIR = os.path.abspath(basepath)

class NtupleToDataframe(object):
    def __init__(self, rootfiles=[]):
        self.rf = rootfiles
        self.dataframe = None

    def add_file(self,rootfile):
        self.rf.append(rootfile)

    def clear_rootfiles(self):
        self.rf = []

    def load_dataframe(self,treename='phaseII'):
        '''Load all values from the given ntuple file into a 
        pandas dataframe.
        
        Args:
            treename (string): name of ROOT tree in files to load into the
            dataframe.  Method will iterate over all rootfiles and form a 
            dataframe with this tree's branches as the data.'''
        dataframe = {}
        for rootf in self.rf:
            f = uproot.open(rootf)
            ftree = f.get(treename)
            branches = ftree.allkeys()
            for branch in branches:
                if branch not in dataframe:
                    dataframe[branch.decode('utf-8')] = list(ftree.get(branch).array())
                elif dataframe[branch.decode('utf-8')][0] is isinstance(np.ndarray):
                    dataframe[branch.decode('utf-8')] = list(dataframe[branch]) + \
                                        list(ftree.get(branch).array())
                else:
                    dataframe[branch.decode('utf-8')] = dataframe[branch] + \
                                        list(ftree.get(branch).array())
                    #dataframe[branch] = np.append(dataframe[branch],
                    #                    ftree.get(branch).array())
        self.dataframe = pd.DataFrame(dataframe)


    def MakeJointPlot(self,xdatalabel,ydatalabel,normed=False):
        if self.dataframe is None:
            print("You must load a dataframe...")
            return
        if isinstance(self.dataframe[xdatalabel][0],list) or \
                isinstance(self.dataframe[ydatalabel][0],list):
            print("You must choose branch types not of jaggedarray to use this method.")
            return
        print("YLABEL ARR IS: " + str(self.dataframe[ydatalabel]))
        g = sns.jointplot(x=xdatalabel,y=ydatalabel,data=self.dataframe,
                kind="hex",stat_func=None).set_axis_labels(xdatalabel,ydatalabel)
        plt.subplots_adjust(left=0.2,right=0.8,
                top=0.95,bottom=0.2)
        g.fig.suptitle("Heatmap of %s and %s"%(xdatalabel,ydatalabel))
        cbar_ax = g.fig.add_axes([0.85,0.2,0.05,0.62])
        plt.colorbar(cax=cbar_ax)
        plt.show()

class AnnieHeatMapMaker(NtupleToDataframe):

    def MakeHitChargeAnglePDF(self,normed=False):
        if self.dataframe is None:
            print("You must load a dataframe...")
            return
        miniframe = {"hitCharges": [], "hitAngles": []}
        print("DATAFRAME SIZE IS: " + str(self.dataframe.size))
        for i in range(len(self.dataframe["trueVtxX"])):
            diType = self.dataframe['digitType'][i]
            PMTs = np.where(diType==0)[0]
            digitQ = self.dataframe['digitQ'][i][PMTs]
            hitangles = self.dataframe['hitAngles'][i][PMTs]
            miniframe["hitCharges"] = miniframe["hitCharges"] + list(digitQ)
            miniframe["hitAngles"] = miniframe["hitAngles"] + list(hitangles)
        miniframe = pd.DataFrame(miniframe)
        ylabel = "Hit angle (degrees)"
        xlabel = "Hit charge"
        g = sns.jointplot(x="hitCharges",y="hitAngles",data=miniframe,
                kind="kde",stat_func=None,color="g").set_axis_labels(xlabel,ylabel)
        #plt.subplots_adjust(left=0.2,right=0.8,
        #        top=0.95,bottom=0.2)
        g.fig.suptitle("Distribution of hit charge and hit angle relative to " +\
                "true muon direction")
        #cbar_ax = g.fig.add_axes([0.85,0.2,0.05,0.62])
        #plt.colorbar(g,cax=cbar_ax)
        plt.show()


    def AddHitAnglesToDataFrame(self):
        if self.dataframe is None:
            print("You must load a dataframe...")
            return
        miniframe = {"hitCharges": [], "hitAngles": []}
        print("DATAFRAME SIZE IS: " + str(self.dataframe.size))
        allhitangles = []
        for i in range(len(self.dataframe["trueVtxX"])):
            digitX = self.dataframe["digitX"][i]
            digitY = self.dataframe["digitY"][i]
            digitZ = self.dataframe["digitZ"][i]
            digitQ = self.dataframe["digitQ"][i]
            trueVtxX = self.dataframe["trueVtxX"][i]
            trueVtxY = self.dataframe["trueVtxY"][i]
            trueVtxZ = self.dataframe["trueVtxZ"][i]
            trueDirX = self.dataframe["trueDirX"][i]
            trueDirY = self.dataframe["trueDirY"][i]
            trueDirZ = self.dataframe["trueDirZ"][i]
            muondir = np.array([trueDirX, trueDirY, trueDirZ])
            hitdirX = digitX-trueVtxX
            hitdirY = digitY-trueVtxY
            hitdirZ = digitZ-trueVtxZ
            hitangles = []
            for j in range(len(hitdirX)):
                hitdir = np.array([digitX[j]-trueVtxX,digitY[j]-trueVtxY,
                               digitZ[j]-trueVtxZ])
                hitdir = hitdir/np.sqrt(np.dot(hitdir,hitdir))
                digitangle = np.arccos(np.dot(hitdir,muondir))
                digitangle = digitangle*180./np.pi
                hitangles.append(digitangle)
            allhitangles.append(np.array(hitangles))
        self.dataframe["hitAngles"] = allhitangles

    def MakePMTHitChargeDataFrame(self):
        if self.dataframe is None:
            print("You must load a dataframe...")
            return
        miniframe = {"hitCharges": [], "hitAngles": []}
        print("DATAFRAME SIZE IS: " + str(self.dataframe.size))
        for i in range(len(self.dataframe["trueVtxX"])):
            diType = self.dataframe['digitType'][i]
            PMTs = np.where(diType==0)[0]
            digitQ = self.dataframe['digitQ'][i][PMTs]
            hitangles = self.dataframe['hitAngles'][i][PMTs]
            miniframe["hitCharges"] = miniframe["hitCharges"] + list(digitQ)
            miniframe["hitAngles"] = miniframe["hitAngles"] + list(hitangles)
        miniframe = pd.DataFrame(miniframe)
        return miniframe

    def MakeYThetaDataFrame(self):
        if self.dataframe is None:
            print("You must load a dataframe...")
            return
        miniframe = {"Y": [], "Theta": []}
        print("DATAFRAME SIZE IS: " + str(self.dataframe.size))
        allhitangles = []
        for i in range(len(self.dataframe["trueVtxX"])):
            digitX = self.dataframe["digitX"][i]
            digitY = self.dataframe["digitY"][i]
            digitZ = self.dataframe["digitZ"][i]
            miniframe["Y"] = miniframe["Y"] + list(digitY)
            miniframe["Theta"] = miniframe["Theta"] +list(hp.XZ_ToTheta(digitX,digitZ))
        return pd.DataFrame(miniframe)
