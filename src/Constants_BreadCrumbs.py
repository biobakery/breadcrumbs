"""
Author: Timothy Tickle
Description: Constants.
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
class Constants_BreadCrumbs():
    """
    Class to hold project constants.
    """

    #References to other projects
    c_strMicropitaProject = "../micropita/src"
    c_strCircladerScript = "../circlader/circlader.py"

    #Character Constants
    COLON = ":"
    COMMA = ","
    ENDLINE = "\n"
    FASTA_ID_LINE_START = ">"
    PATH_SEP = "/"
    QUOTE = "\""
    TAB = '\t'
    WHITE_SPACE = " "
    PIPE = "|"

    #Default values for missing data in the Abundance Table
    c_strEmptyAbundanceData = "0"
    c_strEmptyDataMetadata = "NA"

    lNAs = list(set(["NA","na","Na","nA",c_strEmptyDataMetadata]))
