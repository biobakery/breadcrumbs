#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Make PCoA of an abundance file
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

import argparse
import os
import sys

from AbundanceTable import AbundanceTable
from PCoA import PCoA

ldFilterLevel = [0.001,0.002,0.005]

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "PCoAAbundance.py",
    description = """Creates PCoAs of two different data sets.""" )

#Arguments
#For table
# Example for Spearman ./compareTables.py --id sample --meta study --label study --metric SPEARMAN test.pcl test.pdf
# Example for Bray-Curtis ./compareTables.py --id sample --meta study --label study test.pcl test.pdf
#argp.add_argument("--compareInput", dest="strFileCompare", metavar="FileToCompare", default = None, help="If this is given, this file and the input file are first combined to one abundance table then plotted.")
argp.add_argument("--id", dest="sIDName", metavar= "Sample ID Metadata Name", default="ID", help="Abundance Table ID")
argp.add_argument("--meta", dest="sLastMetadataName", metavar= "Last Metadata Name", help="Last metadata name")
argp.add_argument("--delim", dest= "cFileDelimiter", action= "store", metavar="File Delimiter", default="\t", help="File delimiter, default tab") 
argp.add_argument("--featureDelim", dest="cFeatureNameDelimiter", action= "store", metavar="Feature Name Delimiter", default="|", help="Feature delimiter") 
argp.add_argument("--label", dest="sLabel", metavar= "Label", default=None, help="Label to paint in the PCoA")
argp.add_argument("--doNorm", dest="fDoNormData", action = "store_true", help="Do Normalize Data")
argp.add_argument("--doSum", dest="fDoSumData", action = "store_true", help="Do Sum Data")
argp.add_argument("--metric", dest="strMetric", metavar = "distance", default = PCoA.c_BRAY_CURTIS, help ="Distance metric to use.")
argp.add_argument("--alpha", dest="dAlpha", metavar="alpha", default = 1.0, help="Transparency for points (1.0 is no tranparency, 0.0 is completely transparent..and invisible).")
argp.add_argument("--cladeLevel", dest="iCladeLevel", metavar= "clade level", default=6, type=int, help="Clade level to reduce to (Number 0 = Root).")
argp.add_argument("strFileAbund", metavar = "Abundance file", help ="Input data file")
argp.add_argument("strOutFile", metavar = "Selection Output File", help ="Output file")

args = argp.parse_args( )

#If the compare file is given then combine
#TODO not compete
args.strFileCompare = False
if args.strFileCompare:

    #File pieces
    sInputFileName, sInputFileExt = os.path.splitext(args.strFileAbund)
    sCompareFileBaseName = os.path.splitext(os.path.split(args.strFileCompare)[1])[0]

    #Read in abundance table
    abndInput = AbundanceTable.funcMakeFromFile(args.strFileAbund,
                             cDelimiter = args.cFileDelimiter,
                             sMetadataID = args.sIDName,
                             sLastMetadata = args.sLastMetadataName,
                             cFeatureNameDelimiter= args.cFeatureNameDelimiter)


    #Read in abundance table
    abndCompare = AbundanceTable.funcMakeFromFile(args.strFileCompare,
                             cDelimiter = args.cFileDelimiter,
                             sMetadataID = args.sIDName,
                             sLastMetadata = args.sLastMetadataName,
                             cFeatureNameDelimiter= args.cFeatureNameDelimiter)

    #Sum, normalize and reduce to genus and no otus
    #Would be nice to have a script in here to standardize feature names (with or without the k__)
    abndInput.funcSumClades()
    abndTable.funcNormalize()
    args.fDoSumData = False
    args.fDoNormData = False

    #Combine files
    #TODO Need to write
    
    #Reset input file name
    args.strFileAbund = sInputFileName+"-"+ sCompareFileBaseName+sInputFileExt

    #Write combined file
    abndCombined.funcWriteToFile(args.strFileAbund)

#Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(args.strFileAbund,
                             cDelimiter = args.cFileDelimiter,
                             sMetadataID = args.sIDName,
                             sLastMetadata = args.sLastMetadataName,
                             cFeatureNameDelimiter= args.cFeatureNameDelimiter)

#Get file pieces for input file
sBaseName, sExt = os.path.splitext(args.strOutFile)

#Normalize if need
if args.fDoSumData:
    abndTable.funcSumClades()

#Sum if needed
if args.fDoNormData:
    abndTable.funcNormalize()

#Get the metadata to paint
lsMetadata = abndTable.funcGetMetadata(args.sLabel)

#Get number of samples
iNumberOfSamples = len(abndTable.funcGetSampleNames())

#Get PCoA object and plot
pcoa = PCoA()

#Feature names
lsFeature = abndTable.funcGetFeatureNames()

#Dict to hold thresholding
#Collect on feature average
dictAverageFeatureAbundance = dict([[sFeature,abndTable.funcGetFeatureSumAcrossSamples(sFeature)/iNumberOfSamples] for sFeature in abndTable.funcGetFeatureNames()])

#Reduce to clade and remove otus
if abndTable.funcReduceFeaturesToCladeLevel(args.iCladeLevel):
    abndTable = abndTable.funcGetWithoutOTUs()

    #For each clade in the data, show bray-curtis PCoA
    for dFilterLevel in ldFilterLevel:

        #Filter on filter average
        sFilteredFeatures = [sFeatureName for sFeatureName in dictAverageFeatureAbundance if dictAverageFeatureAbundance.get(sFeatureName,100) <= dFilterLevel]
        abndFeature = abndTable.funcGetFeatureAbundanceTable(lsFeature)

        pcoa.loadData(abndTable,True)
        pcoa.run(tempDistanceMetric=args.strMetric, iDims=2)
        pcoa.plotList(lsLabelList = lsMetadata,
            strOutputFileName = sBaseName+"-c"+str(args.iCladeLevel)+"-f"+str(dFilterLevel)+sExt,
            iSize=20,
            dAlpha=args.dAlpha,
            charForceColor=None,
            charForceShape=None,
            fInvert=False,
            iDim1=1,
            iDim2=2)
else:
    print "Error received when reducing to level, did not occur."
