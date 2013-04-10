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

import csv
import sys
import argparse
import os
from AbundanceTable import AbundanceTable

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "scriptManipulateTable.py",
    description = """Performs common manipulations on tables.\nExample: python scriptManipulateTable.py -i TID -l STSite Test.pcl""" )

#Arguments
#Describe table
argp.add_argument("-i","--id", dest="sIDName", default="ID", help="Abundance Table ID")
argp.add_argument("-l","--meta", dest="sLastMetadataName", help="Last metadata name")
argp.add_argument("-d","--fileDelim", dest= "cFileDelimiter", action= "store", default="\t", help="File delimiter, default tab")
argp.add_argument("-f","--featureDelim", dest= "cFeatureDelimiter", action= "store", default="|", help="Feature (eg. bug or function) delimiter, default '|'")

#Checked x 2
argp.add_argument("-n","--doNorm", dest="fNormalize", action="store_true", default=False, help="Flag to turn on normalization")
argp.add_argument("-s","--doSum", dest="fSum", action="store_true", default=False, help="Flag to turn on summation")

#Unsupervised filtering
argp.add_argument("-A","--doFilterAbundance", dest="strFilterAbundance", action="store", default=None, help="Turns on filtering by abundance (remove features that do not have a minimum abundance in a minimum number of samples); Should be a real number and an integer in the form 'minAbundance,minSamples', (should be performed on a normalized file).")
argp.add_argument("-P","--doFilterPercentile", dest="strFilterPercentile", action="store", default=None, help="Turns on filtering by percentile Should be two numbers between 0 and 1 in the form 'percentile,percentage'. (should be performed on a normalized file).")
argp.add_argument("-O","--doFilterOccurrence", dest="strFilterOccurence", action="store", default=None, help="Turns on filtering by occurrence. Should be two integers in the form 'minSequence,minSample' (should NOT be performed on a normalized file).")
#argp.add_argument("-D","--doFilterDeviation", dest="dCuttOff", action="store", type=float, default=None, help="Flag to turn on filtering by standard deviation (should NOT be performed on a normalized file).")

#Change bug membership
argp.add_argument("-t","--makeTerminal", dest="fMakeTerminal", action="store_true", default=False, help="Works reduces the file to teminal features in the original file.")
argp.add_argument("-u","--reduceOTUs", dest="fRemoveOTUs", action="store_true", default=False, help="Remove otu entries from file.")
argp.add_argument("-c","--reduceToClade", dest="iClade", action="store", type=int, default=None, help="Specify a level of clade to reduce to [].")
argp.add_argument("-b","--reduceToFeatures", dest="strFeatures", action="store", default=None, help="Reduce measurements to certain features (bugs or functions). This can be a comma delimited string (of atleast 2 bugs) or a file.")

#Manipulate based on metadata
#Checked
argp.add_argument("-y","--stratifyBy", dest="strStratifyBy", action="store", default=None, help="Metadata to stratify tables by.")
argp.add_argument("-r","--removeMetadata", dest="strRemoveMetadata", action="store", default=None, help="Remove samples of this metadata and value (format comma delimited string with metadata id first and the values to remove after 'id,lvalue1,value2').")

#Manipulate lineage
#Checked
argp.add_argument("-x","--doPrefixClades", dest="fPrefixClades", action="store_true", default=False, help="Flag to turn on adding prefixes to clades to better identify them, for example s__ will be placed infront of each species.")

#Combine tables
#argp.add_argument("-m","--combineIntersect", dest="fCombineIntersect", action="store_true", default=False, help="Combine two tables including only common features/metadata (intersection).")
#argp.add_argument("-e","--combineUnion", dest="fCombineUnion", action="store_true", default=False, help="Combine two tables (union).")

#Checked
argp.add_argument("-o","--output", dest="strOutFile", action="store", default=None, help="Indicate output pcl file.")
argp.add_argument("strFileAbund", help ="Input data file")


args = argp.parse_args( )

# Creat output file if needed.
if not args.strOutFile:
  args.strOutFile = os.path.splitext(args.strFileAbund)[0]+"-mod.pcl"
lsPieces = os.path.splitext(args.strOutFile)

#List of abundance tables
lsTables = []

#Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(xInputFile=args.strFileAbund,
                                            cDelimiter = args.cFileDelimiter,
                                            sMetadataID = args.sIDName,
                                            sLastMetadata = args.sLastMetadataName,
                                            lOccurenceFilter = None,
                                            cFeatureNameDelimiter=args.cFeatureDelimiter,
                                            xOutputFile = args.strOutFile)

#TODO Check filtering, can not have some filtering together

# Make feature list
lsFeatures = []
if args.strFeatures:
  print "Get features not completed"
#  if "," in args.strFeatures:
#    lsFeatures = args.strFeatures.split(",")
#    print "ManipulateTable::Reading in feature list "+str(len(lsFeatures))+"."
#  else:
#    csvr = csv.reader(open(args.strFeatures, "rU"))
#    print "ManipulateTable::Reading in feature file "+args.strFeatures+"."
#    for lsLine in csvr:
#      lsFeatures.extend(lsLine)

lsTables.append(abndTable)

# Do summing
#Sum if need
if args.fSum:
  for abndTable in lsTables:
    print "ManipulateTable::"+abndTable.funcGetName()+" had "+str(len(abndTable.funcGetFeatureNames()))+" features before summing."
    fResult = abndTable.funcSumClades()
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" was summed."
      print "ManipulateTable::"+abndTable.funcGetName()+" has "+str(len(abndTable.funcGetFeatureNames()))+" features after summing."
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
  lsTerminalTables = []
  for abndTable in lsTables:
    print "ManipulateTable::"+abndTable.funcGetName()+" had "+str(len(abndTable.funcGetFeatureNames()))+" features before making terminal."
    abndTable = abndTable.funcGetFeatureAbundanceTable(abndTable.funcGetTerminalNodes())
    if abndTable:
      print "ManipulateTable::"+abndTable.funcGetName()+" has "+str(len(abndTable.funcGetFeatureNames()))+" terminal features."
      lsTerminalTables.append(abndTable)
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" was not made terminal."
  lsTables = lsTerminalTables

if args.fRemoveOTUs:
  lsNotOTUs = []
  for abndTable in lsTables:
    print "ManipulateTable::"+abndTable.funcGetName()+" had "+str(len(abndTable.funcGetFeatureNames()))+" features before removing OTUs."
    abndTable = abndTable.funcGetWithoutOTUs()
    if abndTable:
      print "ManipulateTable::"+abndTable.funcGetName()+" had OTUs removed and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
      lsNotOTUs.append(abndTable)
    else:
      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" OTUs were not removed."
  lsTables = lsNotOTUs

if args.iClade:
  for abndTable in lsTables:
    fResult = abndTable.funcReduceFeaturesToCladeLevel(args.iClade)
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" was reduced to clade level "+str(args.iClade)+"."
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
      fResult = abndTable.funcFilterAbundanceByPercentile(dPercentileCutOff = float(dPercentile), dPercentageAbovePercentile = float(dPercentage))
      if fResult:
        print "ManipulateTable::"+abndTable.funcGetName()+" has been reduced by percentile and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
      else:
        print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" could not be reduced by percentile."
    else:
      print "ManipulateTable::"+abndTable.funcGetName()+" was NOT normalized and so the percentile filter is invalid, please indicate to normalize the table."

# Filter on abundance (should go after filter on percentile because the filter on percentile 
# needs the full distribution of features in a sample
if args.strFilterAbundance:
  dAbundance,iMinSamples = args.strFilterAbundance.split(",")
  dAbundance = float(dAbundance)
  iMinSamples = int(iMinSamples)
  for abndTable in lsTables:
    if abndTable.funcIsNormalized():
      fResult = abndTable.funcFilterAbundanceByMinValue(dMinAbundance=dAbundance,iMinSamples=iMinSamples)
      if fResult:
        print "ManipulateTable::"+abndTable.funcGetName()+" has been reduced by minimum relative abundance value and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
      else:
        print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" could not be reduced by percentile."
    else:
      print "ManipulateTable::"+abndTable.funcGetName()+" was NOT normalized and so the abundance filter is invalid, please indicate to normalize the table."

#if args.dCuttOff:
#  print "Standard deviation filtering not completed"
#  for abndTable in lsTables:
#    abndTable.funcFilterFeatureBySD(dMinSDCuttOff=args.dCuttOff)
#    if fResult:
#      print "ManipulateTable::"+abndTable.funcGetName()+" has been reduced by standard deviation and now has "+str(len(abndTable.funcGetFeatureNames()))+" features."
#    else:
#      print "ManipulateTable::ERROR. "+abndTable.funcGetName()+" could not be reduced by standard devation."

# Need to normalize again after abundance data filtering given removing features breaks the normalization
# This happends twice because normalization is required to make the abundance data to filter on ;-)
# Normalize if needed
if args.fNormalize:
  for abndTable in lsTables:
    fResult = abndTable.funcNormalize()
    if fResult:
      print "ManipulateTable::"+abndTable.funcGetName()+" was normalized after filtering on abundance data."

#Manipulate lineage
if args.fPrefixClades:
  for abndTable in lsTables:
    fResult = abndTable.funcAddCladePrefixToFeatures()
    if fResult:
      print "ManipulateTable::Clade Prefix was added to "+abndTable.funcGetName()
    else:
      print "ManipulateTable::ERROR. Clade Prefix was NOT added to "+abndTable.funcGetName()

#Manipulate based on metadata
if args.strStratifyBy:
  labndStratifiedTables = []
  for abndTable in lsTables:
    labndResult = abndTable.funcStratifyByMetadata(strMetadata=args.strStratifyBy)
    print "ManipulateTable::"+abndTable.funcGetName()+" was stratified by "+args.strStratifyBy+" in to "+str(len(labndResult))+" tables."
    labndStratifiedTables.extend(labndResult)
  lsTables = labndStratifiedTables

if len(lsTables) == 1:
  lsTables[0].funcWriteToFile(args.strOutFile)
else:
  iIndex = 1
  for abndManTable in lsTables:
    abndManTable.funcWriteToFile(lsPieces[0]+str(iIndex)+lsPieces[1])
    iIndex = iIndex + 1
