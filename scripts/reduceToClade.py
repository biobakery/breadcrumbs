#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Reduce clades of an abundance table to a certain clade level
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
argp = argparse.ArgumentParser( prog = "ReduceToClade.py",
    description = """Reduce clades of an abundance table to a certain clade level.""" )

#Arguments
#For table
argp.add_argument("-id", dest="sIDName", metavar= "Sample ID Metadata Name", default="ID", help="Abundance Table ID")
argp.add_argument("-meta", dest="sLastMetadataName", metavar= "Last Metadata Name", help="Last metadata name")
argp.add_argument("-fDelim", dest= "cFileDelimiter", action= "store", metavar="File Delimiter", default="\t", help="File delimiter, default tab") 
argp.add_argument("-featureDelim", dest="cFeatureNameDelimiter", action= "store", metavar="Feature Name Delimiter", default="|", help="Feature delimiter") 
argp.add_argument("-cladeLevel", dest="iCladeLevel", metavar= "clade level", type=int, help="Clade level to reduce to (Number 0 = Root).")
argp.add_argument("strFileAbund", metavar = "Abundance file", help ="Input data file")
argp.add_argument("strOutFile", metavar = "Selection Output File", help ="Output file")

args = argp.parse_args( )
#Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(args.strFileAbund,
                             cDelimiter = args.cFileDelimiter,
                             sMetadataID = args.sIDName,
                             sLastMetadata = args.sLastMetadataName,
                             cFeatureNameDelimiter= args.cFeatureNameDelimiter)

#Reduce to clade level
if abndTable.funcReduceFeaturesToCladeLevel(args.iCladeLevel):
    abndTable = abndTable.funcGetWithoutOTUs()
    abndTable.funcWriteToFile(xOutputFile=args.strOutFile)
else:
    print "Error received when reducing to level, did not occur."
