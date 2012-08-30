"""
Author: Timothy Tickle
Description: Project constants.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

##
#Used to test the FileIO class
class ConstantsBreadCrumbs():
    """
    Class to hold project constants.
    """

    #Character Constants
    c_strComma = ','
    c_strColon = ':'
    c_strConfigFileHeaderChar = '['
    c_strConfigFileCommentChar = '#'
    c_strEndline = '\n'
    c_strExtDelim = '.'
    c_cFastaIDLineStart = '>'
    c_strPathDelim = '/'
    c_cPipe = '|'
    c_cQuote = '\"'
    c_cTab = '\t'
    c_strWhiteSpace = ' '

    c_strBreadCrumbsSVMSpace = c_strWhiteSpace

    #Default values for missing data in the Abundance Table
    c_strEmptyAbundanceData = "0"
    c_strEmptyDataMetadata = "NA"
    c_strSVMNoSample = "-"

    lNAs = list(set(["NA","na","Na","nA",c_strEmptyDataMetadata]))

    #TODO remove
    #Reference to circlader
    c_strCircladerScript = "circlader/circlader.py"

    #AbundanceTable
    #Suffix given to a file that is check with the checkRawDataFile method
    OUTPUT_SUFFIX = "-checked.pcl"

    def __init__(self):
      pass
