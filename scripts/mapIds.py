#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Maps ids with a mapper file.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

import argparse
import csv
from operator import itemgetter
import re
import sys

#Constants

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "mapIds.py",
    description = """Maps ids in a pcl file.""" )

#Arguments
argp.add_argument("strInputPCL", metavar = "SummaryFile", type = argparse.FileType("r"), help ="Input PCL to be mapped, currently assumes the ids to be in the first row.")
argp.add_argument("strInputMapFile", metavar = "HeaderFile", type = argparse.FileType("r"), help ="Input mapper file, currently assumes currentID\tIDToChangeTo.")
argp.add_argument("strOutputPCL", metavar = "AnnotationFile", type = argparse.FileType("w"), help ="Output file with ids mapped.")

args = argp.parse_args( )

#Read in the summary file and transform to class based descriptions
csvPCL = open(args.strInputPCL,'r') if isinstance(args.strInputPCL, str) else args.strInputPCL
fPCL = csv.reader(csvPCL, delimiter="\t")

#Get header row
lsHeaderRow = fPCL.next()

#Read in mappings
csvMap = open(args.strInputMapFile,'r') if isinstance(args.strInputMapFile, str) else args.strInputMapFile
dictMaps = dict(csv.reader(csvMap, delimiter="\t"))

#Map ids
lsNewIds = [dictMaps.get(sID,sID) for sID in lsHeaderRow]

#Write out ids
#Write out the rest of the file
csvOut = open(args.strOutputPCL,'w') if isinstance(args.strOutputPCL, str) else args.strOutputPCL
fOut = csv.writer(csvOut, delimiter="\t")
fOut.writerows([lsNewIds])
fOut.writerows(fPCL)
csvPCL.close()
csvMap.close()
csvOut.close()
