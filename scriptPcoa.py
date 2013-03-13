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

import sys
import argparse
from AbundanceTable import AbundanceTable
import os
from PCoA import PCoA

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "scriptPcoa.py",
    description = """PCoAs an abundance file given a metadata.\nExample:python scriptPcoa.py -i TID -l STSite""" )

#Arguments
#For table
argp.add_argument("-i","--id", dest="sIDName", default="ID", help="Abundance Table ID")
argp.add_argument("-l","--meta", dest="sLastMetadataName", help="Last metadata name")
argp.add_argument("-d","--fDelim", dest= "cFileDelimiter", action= "store", default="\t", help="File delimiter, default tab")
argp.add_argument("-f","--featureDelim", dest="cFeatureNameDelimiter", action= "store", metavar="Feature Name Delimiter", default="|", help="Feature delimiter") 

argp.add_argument("-n","--doNorm", dest="fDoNormData", action="store_true", default=False, help="Flag to turn on normalization")
argp.add_argument("-s","--doSum", dest="fDoSumData", action="store_true", default=False, help="Flag to turn on summation")

argp.add_argument("-p","--paint", dest="sLabel", metavar= "Label", default=None, help="Label to paint in the PCoA")
argp.add_argument("-m","--metric", dest="strMetric", metavar = "distance", default = PCoA.c_BRAY_CURTIS, help ="Distance metric to use.")
argp.add_argument("-o","--outputFile", dest="strOutFile", metavar= "outputFile", default=None, help="Specify the path for the output figure.")

argp.add_argument("strFileAbund", metavar = "Abundance file", help ="Input data file")


args = argp.parse_args( )

#Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(args.strFileAbund,
                             cDelimiter = args.cFileDelimiter,
                             sMetadataID = args.sIDName,
                             sLastMetadata = args.sLastMetadataName,
                             cFeatureNameDelimiter= args.cFeatureNameDelimiter)

#Normalize if need
if args.fDoSumData:
  abndTable.funcSumClades()

#Sum if needed
if args.fDoNormData:
  abndTable.funcNormalize()

#Get the metadata to paint
lsKeys = abndTable.funcGetMetadataCopy().keys() if not args.sLabel else [args.sLabel]

#Get pieces of output file
if not args.strOutFile:
  args.strOutFile = os.path.splitext(args.strFileAbund)[0]+"-pcoa.pdf"
lsFilePieces = os.path.splitext(args.strOutFile)

lsOutputFiles = [lsFilePieces[0]+"-"+sKey+lsFilePieces[1] for sKey in lsKeys] if not args.sLabel else [args.strOutFile]

for iIndex in xrange(len(lsKeys)):
	lsMetadata = abndTable.funcGetMetadata(lsKeys[iIndex])

	#Get PCoA object and plot
	pcoa = PCoA()
	pcoa.loadData(abndTable,True)
	pcoa.run(tempDistanceMetric=args.strMetric, iDims=2)
	pcoa.plotList(lsLabelList = lsMetadata,
              strOutputFileName = lsOutputFiles[iIndex],
              iSize=20,
              dAlpha=1.0,
              charForceColor=None,
              charForceShape=None,
              fInvert=False,
              iDim1=1,
              iDim2=2)
