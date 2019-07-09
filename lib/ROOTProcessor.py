import uproot
import numpy as np

class ROOTProcessor(object):
    def __init__(self,treename="phaseII"):
        print("ROOTProcessor initializing")
        self.treename = treename
        self.processed_data = {}

    def getProcessedData(self):
        if not self.processed_data:
            print("processed data array empty.  Did you process any data?")
        return self.processed_data

    def clearProcessedData(self):
        '''
        Empties the current processed_data dictionary.
        '''
        self.processed_data = {}

    def setTreeName(self,treename="phaseII"):
        '''
        Set the tree of data to be accessed from loaded ROOT files.

        Input:
            treename [string]
            Tree in ROOT datafiles to be loaded.
        '''
        self.treename = treename

    def addROOTFile(self,rootfile,branches_to_get=None):
        '''
        Opens the root ntuple file at the path given and append it's data
        to the current data dictionary.
        
        Input:
            rootfile [string]
            Path to rootfile to open.

            branches_to_get [array]
            Array of data variables to load into the dictionary.  Use to
            specify specific data types in ROOT file to collect.
        '''
        f = uproot.open(rootfile)
        ftree = f.get(self.treename)
        all_data = ftree.keys()
        print(all_data)
        for dattype in all_data:
            if branches_to_get is not None:
                if dattype.decode('utf-8') not in branches_to_get:
                    continue
            print("LET'S PROCESS SOME ENTRIES")
            thistype_processed = ftree.get(dattype).array()
            self._appendProcessedEntry(dattype,thistype_processed)
            

    def _appendProcessedEntry(self, dattype, thistype_processed):
        '''
        Appends processed data to the current processed_data dictionary.

        Inputs:
            dattype [string]
            key associated with the data location in self.processed_data.

            thistype_processed [array]
            numpy array of data pulled from a ROOT file using uproot.
        '''
        dattype_key = dattype.decode('utf-8')
        print(dattype_key)
        if dattype_key in self.processed_data:
            thistype_proclist = thistype_processed.tolist()
            self.processed_data[dattype_key] = self.processed_data[dattype] + \
                                           thistype_proclist
        else:
            self.processed_data[dattype] = []
            thistype_proclist = thistype_processed.tolist()
            self.processed_data[dattype_key] = self.processed_data[dattype] + \
                                           thistype_proclist
