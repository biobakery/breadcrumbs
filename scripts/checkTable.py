#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Checks tables
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@hsph.harvard.edu"
__status__ = "Development"

import sys
import argparse
from AbundanceTable import AbundanceTable

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "checkTable.py",
    description = """Checks tables.""" )

#Arguments
#For table
argp.add_argument("-meta", dest="sLastMetadataName", metavar= "Last Metadata Name", help="Last metadata name")
argp.add_argument("-fDelim", dest= "cFileDelimiter", action= "store", metavar="File Delimiter", default="\t", help="File delimiter, default tab") 
argp.add_argument("strFileAbund", metavar = "Abundance file", help ="Input data file")
argp.add_argument("strOutFile", metavar = "Selection Output File", help ="Output file")

args = argp.parse_args( )

#Read in abundance table
abndTable = AbundanceTable.funcCheckRawDataFile(strReadDataFileName = args.strFileAbund,
                                                sLastMetadataName = args.sLastMetadataName,
                                                strOutputFileName = args.strOutFile,
                                                cDelimiter = args.cFileDelimiter)
