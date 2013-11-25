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

lsTypeChoices = [ConstantsBreadCrumbs.c_strPCLFile, ConstantsBreadCrumbs.c_strBiomFile]

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "convertBetweenBIOMAndPCL.py",
    description = """Converts a PCL file to a BIOM file and visa versa.""" )

#Arguments
#For table
argp.add_argument("-i","--id", dest = "sID", default = None, metavar = "Sample ID", help = "The metadata indicating the sample ID.")
argp.add_argument("-l","--meta", dest = "sLastMetadataName", default = None, metavar = "Last Metadata Name", help = "The last listed metadata before the first data measurement in the pcl file or to be in the pcl file.")
argp.add_argument("-r","--rowMetadataID", dest = "sLastMetadataRow", default = None,  metavar = "Last Row Metadata Column", help = "If row metadata is present in a PCL file, what is the id of the last row metadata column (most right column that contains row metadata). PCL file only.")
argp.add_argument("-f", "--delim", dest = "cFileDelimiter", action = "store", metavar = "File Delimiter", default = "\t", help = "File delimiter, default tab.") 
argp.add_argument("-p", "--pipe", dest = "fPipe", action = "store_true", help = "If this flag is used, the input is read from standard in and wrote to standard out. When using this option please also specify the format to write in with -f.")
argp.add_argument("-m", "--format", dest = "strFormat", choices = lsTypeChoices, action = "store", metavar = "save_format", default = None, help = "The file format to use when saving, this is not required unless using piping. When using file paths, file format will be assumed. If this is set, it takes priority over automatic guessing.")
argp.add_argument("strFileAbund", metavar = "Abundance file", nargs = "?", help = "Input data file")
argp.add_argument("strOutputFile", metavar = "Selection Output File", nargs = "?", help ="Output file")

args = argp.parse_args( )

# File type to convert to
strOutputFileType = args.strFormat

# Read file format
strInputFileFormat = None
if args.strFormat:
  strInputFileFormat = ConstantsBreadCrumbs.c_strBiomFile if args.strFormat == ConstantsBreadCrumbs.c_strPCLFile else ConstantsBreadCrumbs.c_strPCLFile

# If piping is used, set the handles correctly
if args.fPipe:
  args.strOutputFile = sys.stdout
  args.strFileAbund = sys.stdin
else:

  # Make the output file name (if not given) and get the type of output file name
  # Change the extension from BIOM to pcl
  if strOutputFileType is None:
    lsFilePieces = os.path.splitext(args.strFileAbund)
  strOutputFileType = ConstantsBreadCrumbs.c_strPCLFile if lsFilePieces[-1]=="."+ConstantsBreadCrumbs.c_strBiomFile else ConstantsBreadCrumbs.c_strBiomFile

  if not args.strOutputFile:
    args.strOutputFile = lsFilePieces[0] + "." + strOutputFileType

# Set the last metadata to the id if not given.
if not args.sLastMetadataName:
  args.sLastMetadataName = args.sID

if strOutputFileType is None:
  print( "Please provide a format to store the file as, using -f." )
else:
  # Read in abundance table
  abndTable = AbundanceTable.funcMakeFromFile(args.strFileAbund, cDelimiter=args.cFileDelimiter, sMetadataID=args.sID, sLastMetadataRow = args.sLastMetadataRow, sLastMetadata=args.sLastMetadataName, strFormat = strInputFileFormat)
  if not abndTable:
    print("Could not create an abundance table from the given file and settings.")
  else:
    abndTable.funcWriteToFile(args.strOutputFile, cDelimiter=args.cFileDelimiter, cFileType=strOutputFileType)
