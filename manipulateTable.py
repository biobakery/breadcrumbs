#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Peforms common manipulations on tables
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
import os
from AbundanceTable import AbundanceTable

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "manipulateTable.py",
    description = """Performs comman manipulations on tables.""" )

#Arguments
#Describe table
argp.add_argument("-i","--id", dest="sIDName", default="ID", help="Abundance Table ID")
argp.add_argument("-l","--meta", dest="sLastMetadataName", help="Last metadata name")
argp.add_argument("-d","--fDelim", dest= "cFileDelimiter", action= "store", default="\t", help="File delimiter, default tab")

argp.add_argument("-n","--doNorm", dest="fNormalize", action="store_true", default=False, help="Flag to turn on normalization")
argp.add_argument("-s","--doSum", dest="fSum", action="store_true", default=False, help="Flag to turn on summation")

#Unsupervised filtering
argp.add_argument("-p","--doFilterPercentile", dest="strFilterPercentile", action="store", default=False, help="Flag to turn on filtering by percentile Should be two numbers between 0 and 1 in the form 'percentile,percentage'. (should be performed on a normalized file).")
argp.add_argument("-o","--doFilterOccurrence", dest="strFilterOccurence", action="store", default=None, help="Flag to turn on filtering by occurence. Should be two integers in the form 'minSequence,minSample' (should NOT be performed on a normalized file).")
argp.add_argument("-v","--doFilterDeviation", dest="dCuttOff", action="store", type=float, default=None, help="Flag to turn on filtering by standard deviation (should NOT be performed on a normalized file).")

#Change bug membership
argp.add_argument("-t","--makeTerminal", dest="fMakeTerminal", action="store_true", default=False, help="Works reduces the file to teminal features in the original file.")
argp.add_argument("-u","--reduceOTUs", dest="fRemoveOTUs", action="store_true", default=False, help="Remove otu entries from file.")
argp.add_argument("-c","--reduceToClade", dest="iClade", action="store", type=int, default=None, help="Specify a level of clade to reduce to [].")
argp.add_argument("-b","--reduceToFeatures", dest="strFeatures", action="store", default=None, help="Reduce measurements to certain features (bugs or functions). This can be a comma delimited string or a file.")

#Manipulate based on metadata
argp.add_argument("-y","--stratifyBy", dest="strStratifyBy", action="store", default=None, help="Metadata to stratify tables by.")
argp.add_argument("-r","--removeMetadata", dest="strRemoveMetadata", action="store", default=None, help="Remove samples of this metadata and value (format comma delimited string with metadata id first and the values to remove after 'id,lvalue1,value2').")

#Manipulate lineage
argp.add_argument("-x","--doPrefixClades", dest="fPrefixClades", action="store_true", default=False, help="Flag to turn on adding prefixes to clades to better identify them, for example s__ will be placed infront of each species.")

#Combine tables
argp.add_argument("-m","--combineIntersect", dest="fCombineIntersect", action="store_true", default=False, help="Combine two tables including only common features/metadata (intersection).")
argp.add_argument("-e","--combineUnion", dest="fCombineUnion", action="store_true", default=False, help="Combine two tables (union).")

#Misc
argp.add_argument("-a","--returnAverageSample", dest="fAverage", action="store_true", default=False, help="If selected a synthetic average sample will be returned, this is a last step so if stratification or reduction to features is selected then average samples will be given for those tables.")


argp.add_argument("strFileAbund", help ="Input data file")
argp.add_argument("strOutFile", help ="Output file")

args = argp.parse_args( )

#List of abundance tables
lsTables = []

#Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(xInputFile=args.strFileAbund,
                                            cDelimiter = args.cFileDelimiter,
                                            sMetadataID = args.sIDName,
                                            sLastMetadata = args.sLastMetadataName,
                                            lOccurenceFilter = None,
                                            cFeatureNameDelimiter="_",
                                            xOutputFile = args.strOutFile)

#TODO Check filtering, can not have some filtering together

# Make feature list
lsFeatures = []
if args.strFeatures:
  if "," in strFeatures:
    lsFeatures = strFeature.split(",")
    print "ManipulateTable::Reading in feature list "+len(lsFeatures)+"."
  else:
    csvr = csv.reader(open(args.strFeatures, "rU"))
    print "ManipulateTable::Reading in feature file "+args.strFeatures+"."
    for lsLine in csvr:
      lsFeatures.extend(lsLine)

# Combine tables
#if fCombineIntersect:
#
#if fCombineUnion:
#
lsTables.append(abndTable)

lsPieces = os.path.splitext(args.strOutFile)

#Manipulate lineage
if args.fPrefixClades:
  for abndTable in lsTables:
    fResult = abndTable.funcAddCladePrefixToFeatures()
    if fResult:
      print "ManipulateTable::Clade Prefix was added to "+abndTable.funcGetName()
    else:
      print "ManipulateTable::ERROR. Clade Prefix was NOT added to "+abndTable.funcGetName()

# Do summing and normalization
#Sum if need
if args.fSum:
  for abndTable in lsTables:
    fResult = abndTable.funcSumClades()
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" was summed."
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" was NOT summed."

# Filter on counts
if args.strFilterOccurence:
  iMinimumSequence,iMinimumSample = args.strFilterOccurence.split(",")
  for abndTable in lsTables:
    if abndTable.funcIsNormalized():
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" is normalized and can not be filtered by occurence. This filter needs counts."
    else:
      fResult = abndTable.funcFilterAbundanceBySequenceOccurence(iMinSequence = int(iMinimumSequence), iMinSamples = int(iMinimumSample))
      if fResult:
        print "ManipulateTable::"+abndTable.funcGetName()+" was filtered by occurence and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
      else:
        print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" was NOT filtered by occurence."

# Change bug membership
if args.fMakeTerminal:
  for abndTable in lsTables:
    abndTable.funcGetTerminalNodes()
    print "ManipulateTable::"+abndTable.funcGetName()+" has "+str(len(abndTable.funcGetFeatureNames()))+" terminal features."

if args.fRemoveOTUs:
  for abndTable in lsTables:
    #TODO check if this persists in list
    abndTable = abndTable.funcGetWithoutOTUs()
    if abndTable:
      print "ManipulateTable::"+abndTable.funcGetName()+" had OTUs removed and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" OTUs were not removed."

if args.iClade:
  for abndTable in lsTables:
    fResult = abndTable.funcReduceFeaturesToCladeLevel(iClade)
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" was reduced to clade level "+str(iClade)+"."
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" was NOT reduced in clade levels."

if args.strFeatures:
  for abndTable in lsTables:
    fResult = abndTable.funcGetFeatureAbundanceTable(lsFeatures)
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" has been reduced to given features and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" could not be reduced to the given list."

if args.strRemoveMetadata:
  lsMetadata = args.strRemoveMetadata.split(",")
  for abndTable in lsTables:
    fResult = abndTable.funcRemoveSamplesByMetadata(sMetadata=lsMetadata[0], lValuesToRemove=lsMetadata[1:])
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" has had samples removed and now has "+str(len(abndTable.funcGetSampleNames()))+" samples."
    else:
      print "ManipulateTable::ERROR. Could not remove samples from "+abndTable.funcGetName()+"."

# Normalize if needed
if args.fNormalize:
  for abndTable in lsTables:
    #TODO
    fResult = abndTable.funcNormalize()
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" was normalized."
    else:
      print "ManipulateTable::"+abndTable.funcGetName()+" was NOT normalized."

# Filter on percentile
if args.strFilterPercentile:
  dPercentile,dPercentage = args.strFilterPercentile.split(",")
  for abndTable in lsTables:
    if abndTable.funcIsNormalized():
      abndTable.funcFilterAbundanceByPercentile(dPercentileCutOff = str(dPercentile), dPercentageAbovePercentile = str(dPercentage))
      if fResult:
        print "ManipulateTable::"+abndTable.funcGetName()+" has been reduced by percentile and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
      else:
        print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" could not be reduced by percentile."
    else:
      print "ManipulateTable::"+abndTable.funcGetName()+" was NOT normalized and so the percentile filter is invalid, please indicate to normalize the table."

if args.dCuttOff:
  for abndTable in lsTables:
    abndTable.funcFilterFeatureBySD(dMinSDCuttOff=args.dCuttOff)
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" has been reduced by standard deviation and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" could not be reduced by standard devation."

#Manipulate based on metadata
if args.strStratifyBy:
  labndStratifiedTables = []
  for abndTable in lsTables:
    labndResult = abndTable.funcStratifyByMetadata(strMetadata=args.strStratifyBy)
    print "ManipulateTable::"+abndTable.funcGetName()+" was stratified by "+args.strStratifyBy+" in to "+str(len(labndResult))+" tables."
    labndStratifiedTables.extend(labndResult)
  lsTables = labndStratifiedTables

# Misc
#TODO
#if fAverage:
#  for abndTable in lsTables:
#    abndTable.funcGetAverageSample()

if len(lsTables) == 1:
  lsTables[0].funcWriteToFile(args.strOutFile)
else:
  iIndex = 1
  for abndManTable in lsTables:
    abndManTable.funcWriteToFile(lsPieces[0]+str(iIndex)+lsPieces[1])
    iIndex = iIndex + 1
