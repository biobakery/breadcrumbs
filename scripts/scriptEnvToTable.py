#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Convert Env file to table
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
import csv


#Set up arguments reader
argp = argparse.ArgumentParser( prog = "scriptEnvToTable.py",
    description = """Convert Env file to table""" )

#Arguments
#For table
argp.add_argument("strEnvFile", metavar = "EnvFile", help ="EnvFile data file")
argp.add_argument("strOutputFile", metavar = "OutputFile", help ="Output File")
args = argp.parse_args( )

hndlReader = csv.reader(open(args.strEnvFile,'rU'), delimiter="\t")

lsListOfIDs = []
lsListOfFeatures = []
dictValues = {}
for lsLine in hndlReader:
  print(lsLine)
  lsListOfIDs.append(lsLine[1])
  lsListOfFeatures.append(lsLine[0])
  tpleKey = tuple([lsLine[1],lsLine[0]])
  if tpleKey in dictValues:
    print("Error:: Duplicate key entries found")
    exit(1)
  dictValues[tpleKey] = lsLine[2]

lsListOfIDs = list(set(lsListOfIDs))
lsListOfFeatures = list(set(lsListOfFeatures))
print(lsListOfIDs)
print(lsListOfFeatures)
hndlWrite = csv.writer(open(args.strOutputFile,'w'), delimiter="\t")
hndlWrite.writerow(["ID"]+lsListOfIDs)
for sFeature in lsListOfFeatures:
  lsFeatureLine = [sFeature]
  for sSample in lsListOfIDs:
    lsFeatureLine.append(dictValues.get(tuple([sSample,sFeature]),0))
  hndlWrite.writerow(lsFeatureLine)
