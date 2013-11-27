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
from breadcrumbs.src.AbundanceTable import AbundanceTable
from breadcrumbs.src.Metric import Metric
import csv
import os
from breadcrumbs.src.PCoA import PCoA

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
argp.add_argument("-m","--metric", dest="strMetric", metavar = "distance", default = "braycurtis", help ="Distance metric to use. Pick from braycurtis, canberra, chebyshev, cityblock, correlation, cosine, euclidean, hamming, spearman, sqeuclidean, unifrac_unweighted, unifrac_weighted")
argp.add_argument("-o","--outputFile", dest="strOutFile", metavar= "outputFile", default=None, help="Specify the path for the output figure.")
argp.add_argument("-D","--DistanceMatrix", dest="strFileDistanceMatrix", metavar= "strFileDistanceMatrix", default=None, help="Specify the path for outputing the distance matrix (if interested). Default this will not output.")
argp.add_argument("-C","--CoordinatesMatrix", dest="strFileCoordinatesMatrix", metavar= "strFileCoordinatesMatrix", default=None, help="Specify the path for outputing the x,y coordinates matrix (Dim 1 and 2). Default this will not output.")

# Unifrac arguments
argp.add_argument("-t","--unifracTree", dest="istrmTree", metavar="UnifracTreeFile", default=None, help="Optional file only needed for UniFrac calculations.")
argp.add_argument("-e","--unifracEnv", dest="istrmEnvr", metavar="UnifracEnvFile", default=None, help="Optional file only needed for UniFrac calculations.")
argp.add_argument("-c","--unifracColor", dest="fileUnifracColor", metavar="UnifracColorFile", default = None, help="A text file indicating the groupings of metadata to color. Each line in the file is a group to color. An example file line would be  'GroupName:ID,ID,ID,ID'")

argp.add_argument("strFileAbund", metavar = "Abundance file", nargs="?", help ="Input data file")

args = argp.parse_args( )

#Read in abundance table
abndTable = None
if args.strFileAbund:
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
lsKeys = None
if abndTable:
  lsKeys = abndTable.funcGetMetadataCopy().keys() if not args.sLabel else [args.sLabel]

#Get pieces of output file
if not args.strOutFile:
  if not args.strFileAbund:
    args.strOutFile = os.path.splitext(os.path.basename(args.istrmEnvr))[0]+"-pcoa.pdf"
  else:
    args.strOutFile = os.path.splitext(os.path.basename(args.strFileAbund))[0]+"-pcoa.pdf"
lsFilePieces = os.path.splitext(args.strOutFile)

# Make PCoA object
# Get PCoA object and plot
pcoa = PCoA()
if(not args.strMetric in [Metric.c_strUnifracUnweighted,Metric.c_strUnifracWeighted]) and abndTable:
  pcoa.loadData(abndTable,True)
# Optional args.strFileDistanceMatrix if not none will force a printing of the distance measures to the path in args.strFileDistanceMatrix
pcoa.run(tempDistanceMetric=args.strMetric, iDims=2, strDistanceMatrixFile=args.strFileDistanceMatrix, istrmTree=args.istrmTree, istrmEnvr=args.istrmEnvr)

# Write dim 1 and 2 coordinates to file
if args.strFileCoordinatesMatrix:
  lsIds = pcoa.funcGetIDs()
  mtrxCoordinates = pcoa.funcGetCoordinates()
  csvrCoordinates = csv.writer(open(args.strFileCoordinatesMatrix, 'w'))
  csvrCoordinates.writerow(["ID","Dimension_1","Dimension_2"])
  for x in xrange(mtrxCoordinates.shape[0]):
    strId = lsIds[x] if lsIds else ""
    csvrCoordinates.writerow([strId]+mtrxCoordinates[x].tolist())

# Paint metadata
if lsKeys:
  for iIndex in xrange(len(lsKeys)):
    lsMetadata = abndTable.funcGetMetadata(lsKeys[iIndex])

    pcoa.plotList(lsLabelList = lsMetadata,
      strOutputFileName = lsFilePieces[0]+"-"+lsKeys[iIndex]+lsFilePieces[1],
      iSize=20,
      dAlpha=1.0,
      charForceColor=None,
      charForceShape=None,
      fInvert=False,
      iDim1=1,
      iDim2=2)

if args.strMetric in [Metric.c_strUnifracUnweighted,Metric.c_strUnifracWeighted]:

  c_sNotGiven = "Not_specified"

  lsIds = pcoa.funcGetIDs()
  lsGroupLabels = [c_sNotGiven for s in lsIds]

  if args.fileUnifracColor:

    # Read color file and make a dictionary to convert ids
    lsColorLines = csv.reader(open(args.fileUnifracColor))
    dictConvertIDToGroup = {}
    for lsLine in lsColorLines:
      if lsLine:
        sGroupID, sFirstID = lsLine[0].split(":")
        dictConvertIDToGroup.update(dict([(sID,sGroupID) for sID in [sFirstID]+lsLine[1:]]))

    lsGroupLabels = [dictConvertIDToGroup.get(sID,c_sNotGiven) for sID in lsIds]

  pcoa.plotList(lsLabelList = lsGroupLabels,
      strOutputFileName = lsFilePieces[0]+"-"+args.strMetric+lsFilePieces[1],
      iSize=20,
      dAlpha=1.0,
      charForceColor=None,
      charForceShape=None,
      fInvert=False,
      iDim1=1,
      iDim2=2)
