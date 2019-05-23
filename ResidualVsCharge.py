import uproot
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize as scp
import pandas as pd
import numpy as np
import sys

def theplot(data):
    sns.set(font_scale=2)
    sns.regplot(data["Residuals"],data["Charges"],fit_reg=False)
    plt.xlabel(r'Hit Residual (ns)',fontsize=24)
    plt.ylabel(r'Charge',fontsize=24)
    plt.title("Charge vs. hit residual for simulated CCQE events \n" +
            r"residual = $t_{digit} - \frac{L}{c} - t_{muon}$",fontsize=30)
    plt.legend(fontsize=22)
    plt.ion()
    plt.show()

if __name__=='__main__':
    f = uproot.open("./RPTest_100.root")
    ftree = f.get("phaseII")
    ftree.items()
    digitT = ftree.get("digitT")
    digitQ = ftree.get("digitQ")
    TrueVtxTime = ftree.get("trueVtxTime")
    TrueVtxX = ftree.get("trueVtxX")
    TrueVtxY = ftree.get("trueVtxY")
    TrueVtxZ = ftree.get("trueVtxZ")
    digitX = ftree.get("digitX")
    digitY = ftree.get("digitY")
    digitType = ftree.get("digitType")
    digitZ = ftree.get("digitZ")
    evn = ftree.get("eventNumber")
    diQ = digitQ.array()
    evnums = evn.array()
    diT = digitT.array()
    diX = digitX.array()
    diY = digitY.array()
    diZ = digitZ.array()
    diType = digitType.array()
    truX = TrueVtxX.array()
    truY = TrueVtxY.array()
    truZ = TrueVtxZ.array()
    trueT = TrueVtxTime.array()
    tot_entries = len(evn)
    PMT_charges = np.array([])
    PMT_SOL_residuals = np.array([])
    PMT_SOL_wat_residuals = np.array([])
    LAPPD_charges = np.array([]) 
    LAPPD_SOL_residuals = np.array([])
    LAPPD_SOL_wat_residuals = np.array([])

    for e in xrange(tot_entries):
        #Want to calculate all hit residuals based on if
        #Light is moving at C or at C/1.33
        typedata = diType[e]
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
        #Append these charges to the total charge array
        PMT_charges = np.append(this_diQ[pmtind],PMT_charges)
        #Calculate speed of light residuals
        L = np.sqrt((this_diX-this_TX)**2 + (this_diY - this_TY)**2 + (this_diZ - this_TZ)**2)
        #Residuals assuming speed of light in vacuum (cm/ns)
        res_SOL = this_diT - (L*1.0/29.97) - this_TT
        pmt_res_SOL = res_SOL[pmtind]
        PMT_SOL_residuals = np.append(pmt_res_SOL,PMT_SOL_residuals)
        lappd_res_SOL = res_SOL[lappdind]
        LAPPD_SOL_residuals = np.append(lappd_res_SOL,LAPPD_SOL_residuals)
        
        #Residuals assuming speed of light in water (cm/ns)
        res_SOL_wat = this_diT - (L*1.42/29.97) - this_TT
        pmt_res_SOL_wat = res_SOL_wat[pmtind]
        PMT_SOL_wat_residuals = np.append(pmt_res_SOL_wat,PMT_SOL_wat_residuals)
        lappd_res_SOL_wat = res_SOL_wat[lappdind]
        LAPPD_SOL_wat_residuals = np.append(lappd_res_SOL_wat,LAPPD_SOL_wat_residuals)
    #now, we produce datasets giving hit residuals vs. charge for PMTs
    thedata = {"Residuals":PMT_SOL_residuals, "Charges":PMT_charges}
    print(len(thedata["Residuals"]))
    print(len(thedata["Charges"]))
    dataPD = pd.DataFrame(thedata)
    theplot(dataPD)
