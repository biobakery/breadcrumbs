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

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "manipulateTable.py",
    description = """Performs comman manipulations on tables.""" )

#Arguments
#Describe table
argp.add_argument("--id", dest="sIDName", default="ID", help="Abundance Table ID")
argp.add_argument("--meta", dest="sLastMetadataName", help="Last metadata name")
argp.add_argument("--fDelim", dest= "cFileDelimiter", action= "store", default="\t", help="File delimiter, default tab")

argp.add_argument("--doNorm", dest="fNormalize", action="store_true", default=False, help="Flag to turn on normalization")
argp.add_argument("--doSum", dest="fSum", action="store_true", default=False, help="Flag to turn on summation")

#Unsupervised filtering
argp.add_argument("--doFilterPercentile", dest="fFilterPercentile", action="store_true", default=False, help="Flag to turn on filtering by percentile (should be performed on a normalized file).")
argp.add_argument("--doFilterOccurrence", dest="fFilterOccurence", action="store_true", default=False, help="Flag to turn on summation (should NOT be performed on a normalized file).")
argp.add_argument("--doFilterDeviation", dest="fFilterDeviation", action="store_true", default=False, help="Flag to turn on summation (should NOT be performed on a normalized file).")

#Change bug membership
argp.add_argument("--makeTerminal", dest="fMakeTerminal", action="store_true", default=False, help="Works reduces the file to teminal features in the original file.")
argp.add_argument("--reduceOTUs", dest="fRemoveOTUs", action="store", default=None, help="Remove otu entries from file.")
argp.add_argument("--reduceToClade", dest="iClade", action="store", default=None, help="Specify a level of clade to reduce to [].")

#Manipulate based on metadata
argp.add_argument("--doPrefixClades", dest="fPrefixClades", action="store_true", default=False, help="Flag to turn on adding prefixes to clades to better identify them, for example s__ will be placed infront of each species.")



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
#Sum if need
if args.fSum:
  abndTable.funcSumClades()

#Normalize if needed
if args.fNormalize:
  abndTable.funcNormalize()

abndTable.funcWriteToFile(args.strOutFile)
