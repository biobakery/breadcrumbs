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
from AbundanceTable import AbundanceTable
from BoxPlot import BoxPlot

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "boxPlot.py",
    description = """Make a box plot from an abundance table.""" )

#Arguments
argp.add_argument("-t","--title", dest="strTitle", action="store", default=None, help="The Text to be the title.")
argp.add_argument("-x","--xaxis", dest="strX", action="store", default=None, help="Text for the x-axis.")
argp.add_argument("-y","--yaxis", dest="strY", action="store", default=None, help="Text for the y-axis.")
argp.add_argument("-c","--color", dest="strColor", action="store", default=None, help="Hex color (including the #).")

argp.add_argument("-r","--invertcolor", dest="fColor", action="store_true", default=False, help="Flag to invert the background to black.")
argp.add_argument("-s","--invertyaxis", dest="fAxis", action="store_true", default=False, help="Flag to invert the y axis.")

argp.add_argument("strY", help="Id of the metadata to be used as the values (Y axis).")
argp.add_argument("strX", help="Id of the metadata to be used as the groups (X axis).")
argp.add_argument("strFileAbund", help ="Input data file")
argp.add_argument("strOutFile", help ="Output file")

args = argp.parse_args( )


#Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(xInputFile=args.strFileAbund,
                                            cDelimiter = args.cFileDelimiter,
                                            sMetadataID = args.sIDName,
                                            sLastMetadata = args.sLastMetadataName,
                                            lOccurenceFilter = None,
                                            cFeatureNameDelimiter="_",
                                            xOutputFile = args.strOutFile)

strTitleParam = args.strTitle if args.strTitle else "Title"
strXAxisParam = args.strX if args.strX else "X Axis"
strYAxisParam = args.strY if args.strY else "Y Axis"
strColorParam = args.strColor if args.strTitle else "#83C8F9"

ly = abndTable.funcGetFeature(strYAxisParam)
lsLabels = abndTable.funcGetMetadata(strXAxisParam)

BoxPlot.funcPlot(ly=ly, lsLabels=lsLabels, strOutputFigurePath=args.strOutFile, strTitle=strTitleParam, strXTitle=strXAxisParam, strYTitle=strYAxisParam, strColor=strColorParam, fJitter=False, fInvert=args.fColor, fInvertY=args.fAxis)

