#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Plots feaures
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
from breadcrumbs.src import BoxPlot
from breadcrumbs.src import Histogram
from breadcrumbs.src import ScatterPlot

def funcPlotBoxPlot(lxVariable1,lxVariable2,fOneIsNumeric):

  ly,lsLabels = [lxVariable1,lxVariable2] if fOneIsNumeric else [lxVariable2,lxVariable1]

  # Group data
  dictGroups = {}
  for iIndex in xrange(len(ly)):
    lsList = dictGroups.get(lsLabels[iIndex],[])
    lsList.append(ly[iIndex])
    dictGroups.setdefault(lsLabels[iIndex],lsList)
  ly = [dictGroups[sKey] for sKey in dictGroups.keys()]
  lsLabels = dictGroups.keys()

  BoxPlot.funcPlot(ly=ly, lsLabels=lsLabels, strOutputFigurePath=args.strOutputFile, strTitle=args.strTitle, strXTitle=args.strX, strYTitle=args.strY, strColor=args.strColor, fJitter=True, fInvert=args.fColor, fInvertY=args.fAxis)


#Set up arguments reader
argp = argparse.ArgumentParser( prog = "scriptBoxPlot.py\nExample: python scriptBoxPlot.py Input.pcl valuesID groupID",
    description = "Make a box plot from an abundance table.")

#Sepecify output if needed
argp.add_argument("-o","--output", dest="strOutputFile", action="store", default=None, help="Output file name.")

# Text annotation
argp.add_argument("-t","--title", dest="strTitle", action="store", default=None, help="Test for the title.")
argp.add_argument("-x","--xaxis", dest="strX", action="store", default=None, help="Text for the x-axis.")
argp.add_argument("-y","--yaxis", dest="strY", action="store", default=None, help="Text for the y-axis.")

# Color options
argp.add_argument("-c","--color", dest="strColor", action="store", default="#83C8F9", help="Fill color as a Hex number (including the #).")
argp.add_argument("-r","--invertcolor", dest="fColor", action="store_true", default=False, help="Flag to invert the background to black.")

# Axis adjustments
argp.add_argument("-s","--invertyaxis", dest="fAxis", action="store_true", default=False, help="Flag to invert the y axis.")

# Required
argp.add_argument("strFileAbund", help ="Input data file")
argp.add_argument("strFeatures", nargs = "+", help="Features to plot (from one to two metadata).")

args = argp.parse_args( )

#Holds the data
lxVariable1 = None
lxVariable2 = None
fOneIsNumeric = False
fTwoIsNumeric = False

strFeatureOneID = args.strFeatures[0]
strFeatureTwoID = None if len(args.strFeatures)<2 else args.strFeatures[1]

# If the output file is not specified, make it up
if not args.strOutputFile:
  lsPieces = os.path.splitext(args.strFileAbund)
  args.strOutputFile = [lsPieces[0],strFeatureOneID]
  if strFeatureTwoID:
    args.strOutputFile = args.strOutputFile+[strFeatureTwoID]
  args.strOutputFile = "-".join(args.strOutputFile+["plotfeature.pdf"])

if not args.strTitle:
  args.strTitle = [strFeatureOneID]
  if strFeatureTwoID:
    args.strTitle = args.strTitle+[strFeatureTwoID]
  args.strTitle = " vs ".join(args.strTitle)

csvReader = csv.reader(open(args.strFileAbund, 'rU') if isinstance(args.strFileAbund,str) else args.strFileAbund, delimiter="\t")

if args.strX is None:
  args.strX = strFeatureOneID

if args.strY is None:
  args.strY = strFeatureTwoID

# Get values and groupings
for lsLine in csvReader:
  if lsLine[0] == strFeatureOneID:
    lxVariable1 = lsLine[1:]
  if not strFeatureTwoID is None:
    if lsLine[0] == strFeatureTwoID:
      lxVariable2 = lsLine[1:]

# Remove NAs
liNAs = [i for i,x in enumerate(lxVariable1) if x.lower() == "na"]
liNAs = set([i for i,x in enumerate(lxVariable1) if x.lower() == "na"]+liNAs)
lxVariable1 = [x for i,x in enumerate(lxVariable1) if not i in liNAs]

if not lxVariable2 is None:
  lxVariable2 = [x for i,x in enumerate(lxVariable2) if not i in liNAs]

# Type variables
if not lxVariable1 is None:
  try:
    float(lxVariable1[0])
    lxVariable1 = [float(xItem) for xItem in lxVariable1]
    fOneIsNumeric = True
  except ValueError:
    pass

if not lxVariable2 is None:
  try:
    float(lxVariable2[0])
    lxVariable2 = [float(xItem) for xItem in lxVariable2]
    fTwoIsNumeric = True
  except ValueError:
    pass

if lxVariable1 is None:
  print("scriptPlotFeature:: Sorry, could not find the feature "+ strFeatureOneID +" in the file "+args.strFileAbund+" .")
elif( (lxVariable2 is None) and (not strFeatureTwoID is None) ):
  print("scriptPlotFeature:: Sorry, could not find the feature "+ strFeatureTwoID +" in the file "+args.strFileAbund+" .")
else:
  # Plot as needed
  if((not lxVariable1 is None ) and (not lxVariable2 is None)):
    if(sum([fOneIsNumeric, fTwoIsNumeric])==0):
      print "scriptPlotFeature:: Error, If plotting 2 variables, atleast 1 should be numeric."
    elif(sum([fOneIsNumeric, fTwoIsNumeric])==1):
      funcPlotBoxPlot(lxVariable1,lxVariable2, fOneIsNumeric=fOneIsNumeric)
    elif(sum([fOneIsNumeric, fTwoIsNumeric])==2):
      ScatterPlot.funcPlot(lxVariable1, lxVariable2, args.strOutputFile, strTitle=args.strTitle, strXTitle=args.strX, strYTitle=args.strY, strColor=args.strColor, fInvert=args.fColor)
  elif(not lxVariable1 is None ):
    if fOneIsNumeric:
      Histogram.funcPlot(lxVariable1, args.strOutputFile, strTitle=args.strTitle, strXTitle=args.strX, strYTitle="Frequency", strColor=args.strColor, fInvert=args.fColor)
    else:
      print "Sorry currently histograms are support for only numeric data."
