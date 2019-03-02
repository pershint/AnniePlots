import ROOT
from ROOT import gROOT
import time

def ConvertMeshgrid(xx,yy,zz,xlabel,ylabel,zlabel,title):
    #Get the range of the histograms (They should both be the same)
    Nbinx = len(xx)
    Nbiny = len(yy[0])
    Xmin = xx.min()
    Xmax = xx.max()
    Ymin = yy.min()
    Ymax = yy.max()
    #Now, we divide one by the other
    theHist = ROOT.TH2D('ConvertedHistogram',title,
            Nbinx,Xmin,Xmax,Nbiny,Ymin,Ymax)
    for i,xrow in enumerate(xx):
        for j,xbin in enumerate(xrow):
            print("XBIN,YBIN: %i,%i"%(i,j))
            print("ZVALUE: %f"%(zz[i][j]))
            theHist.SetBinContent(i,j,zz[i][j])
    theHist.SetXTitle(xlabel)
    theHist.SetYTitle(ylabel)
    theHist.SetZTitle(zlabel)
    return theHist
    
    #Eventually we could save these, but nah for now
