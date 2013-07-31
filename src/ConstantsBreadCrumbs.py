"""
Author: Timothy Tickle
Description: Project constants.
"""

#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#####################################################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "MIT"
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
    c_matrixFileDelim = '\t'

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

	######################################################################
	# Constants related to biom import and export files                  #
	######################################################################
    c_biom = "biom"							#Suffix of biom files
    c_BiomTaxData = "BiomTaxData"
    c_Metadata = "Metadata"
    c_metadata_lowercase = "metadata"
    c_sLastMetadata = "sLastMetadata"
    c_columns = "columns"	
    c_ascii = "ascii"	
    c_ignore = "ignore"	
    c_Dtype = "Dtype"	
    c_ID = "ID"	
    c_id_lowercase = "id"	
    c_f4 = "f4"		
    c_biom_file_generated_by = "AbundanceTable conversion program"
    c_strPCLFile = "pcl"
    c_strBiomeFile = "biom"

    def __init__(self):
      pass
