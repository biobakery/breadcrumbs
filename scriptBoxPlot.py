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
import csv
import os
from BoxPlot import BoxPlot

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "scriptBoxPlot.py\nExample: python scriptBoxPlot.py Input.pcl valuesID groupID",
    description = "Make a box plot from an abundance table.")

#Sepecify output if needed
argp.add_argument("-o","--output", dest="strOutputFile", action="store", default=None, help="Output file name.")

# Text annotation
argp.add_argument("-t","--title", dest="strTitle", action="store", default="Title", help="Test for the title.")
argp.add_argument("-x","--xaxis", dest="strX", action="store", default="X axis", help="Text for the x-axis.")
argp.add_argument("-y","--yaxis", dest="strY", action="store", default="Y axis", help="Text for the y-axis.")

# Color options
argp.add_argument("-c","--color", dest="strColor", action="store", default="#83C8F9", help="Fill color as a Hex number (including the #).")
argp.add_argument("-r","--invertcolor", dest="fColor", action="store_true", default=False, help="Flag to invert the background to black.")

# Axis adjustments
argp.add_argument("-s","--invertyaxis", dest="fAxis", action="store_true", default=False, help="Flag to invert the y axis.")

# Required
argp.add_argument("strFileAbund", help ="Input data file")
argp.add_argument("strValuesID", help="Id of the metadata to be used as the values (Y axis).")
argp.add_argument("strGroupsID", help="Id of the metadata to be used as the groups (X axis).")

args = argp.parse_args( )

# If the output file is not specified, make it up
if not args.strOutputFile:
  lsPieces = os.path.splitext(args.strFileAbund)
  args.strOutputFile = lsPieces[0]+"-boxplot.pdf"

csvReader = csv.reader(open(args.strFileAbund, 'rU') if isinstance(args.strFileAbund,str) else args.strFileAbund, delimiter="\t")
ly = None
lsLabels = None

# Get values and groupings
for lsLine in csvReader:
  print lsLine[0]
  if lsLine[0] == args.strValuesID:
    ly = [float(strValue) for strValue in lsLine[1:]]
  if lsLine[0] == args.strGroupsID:
    lsLabels = lsLine[1:]

# Group data
dictGroups = {}
for iIndex in xrange(len(ly)):
  lsList = dictGroups.get(lsLabels[iIndex],[])
  lsList.append(ly[iIndex])
  dictGroups.setdefault(lsLabels[iIndex],lsList)
ly = [dictGroups[sKey] for sKey in dictGroups.keys()]
lsLabels = dictGroups.keys()

BoxPlot.funcPlot(ly=ly, lsLabels=lsLabels, strOutputFigurePath=args.strOutputFile, strTitle=args.strTitle, strXTitle=args.strX, strYTitle=args.strY, strColor=args.strColor, fJitter=False, fInvert=args.fColor, fInvertY=args.fAxis)

