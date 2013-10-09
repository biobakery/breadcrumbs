#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Converts between BIOM and PCL files. If a PCL file is read, an equivalent BIOM file will be written; if a BIOM file is read, an equivalent pcl file will be written.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2013"
__credits__ = ["Timothy Tickle","George Weingart"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@hsph.harvard.edu"
__status__ = "Development"

from AbundanceTable import AbundanceTable
import argparse
from ConstantsBreadCrumbs import ConstantsBreadCrumbs
import os
import sys


#Set up arguments reader
argp = argparse.ArgumentParser( prog = "convertBetweenBIOMAndPCL.py",
    description = """Converts a PCL file to a BIOM file and visa versa.""" )

#Arguments
#For table
argp.add_argument("-i","--id", dest="sID", default = None, metavar= "Sample ID", help="The metadata indicating the sample ID.")
argp.add_argument("-l","--meta", dest = "sLastMetadataName", default = None, metavar= "Last Metadata Name", help="The last listed metadata before the first data measurement in the pcl file or to be in the pcl file.")
argp.add_argument("-r","--rowMetadataID", dest = "sLastMetadataRow", default = None,  metavar = "Last Row Metadata Column", help = "If row metadata is present in a PCL file, what is the id of the last row metadata column (most right column that contains row metadata). PCL file only.")
argp.add_argument("-f","--delim", dest = "cFileDelimiter", action= "store", metavar="File Delimiter", default="\t", help="File delimiter, default tab") 
argp.add_argument("strFileAbund", metavar = "Abundance file", help ="Input data file")
argp.add_argument("strOutputFile", default = "", nargs="?", metavar = "Selection Output File", help ="Output file")

args = argp.parse_args( )

# Make the output file name (if not given) and get the type of output file name
# Change the extension from BIOM to pcl
lsFilePieces = os.path.splitext(args.strFileAbund)
strOutputFileType = ConstantsBreadCrumbs.c_strPCLFile if lsFilePieces[-1]=="."+ConstantsBreadCrumbs.c_strBiomFile else ConstantsBreadCrumbs.c_strBiomFile

if not args.strOutputFile:
  args.strOutputFile = lsFilePieces[0] + "." + strOutputFileType

# Set the last metadata to the id if not given.
if not args.sLastMetadataName:
  args.sLastMetadataName = args.sID

# Read in abundance table
abndTable = AbundanceTable.funcMakeFromFile(args.strFileAbund, cDelimiter=args.cFileDelimiter, sMetadataID=args.sID, sLastMetadataRow = args.sLastMetadataRow, sLastMetadata=args.sLastMetadataName, xOutputFile=args.strOutputFile)
if not abndTable:
  print("Could not create an abundance table from the given file and settings.")
else:
  abndTable.funcWriteToFile(args.strOutputFile, cDelimiter=args.cFileDelimiter, cFileType=strOutputFileType)
