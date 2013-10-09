"""
Author: Timothy Tickle
Description: Base class for ordination plots.
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
__copyright__ = "Copyright 2013"
__credits__ = ["Timothy Tickle"]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
import AbundanceTable
from ConstantsFiguresBreadCrumbs import ConstantsFiguresBreadCrumbs
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from UtilityMath import UtilityMath
from ValidateData import ValidateData

class Ordination:
    """
    Base class for ordination methods and plots.
    """

    def __init__(self):
      # Rows = Samples
      self.dataMatrix = None
      self.isRawData = None
      self.lsIDs = []
      self.dataProcessed = None

    #Happy path tested
    def loadData(self, xData, fIsRawData):
        """
        Loads data into the object (given a matrix or an abundance table)
        Data can be the Abundance Table to be converted to a distance matrix or a distance matrix
        If it is the AbundanceTable, indicate that it is rawData (tempIsRawData=True)
        If it is the distance matrix already generated indicate (tempIsRawData=False)
        and no conversion will occur in subsequent methods.

        :params xData: AbundanceTable or Distance matrix . Taxa (columns) by samples (rows)(lists)
        :type: AbundanceTable or DistanceMatrix
        :param fIsRawData: Indicates if the xData is an AbudanceTable (True) or distance matrix (False; numpy array)
        :type: boolean
        :return boolean: indicator of success (True=Was able to load data)
        """

        if fIsRawData:
            #Read in the file data to a numpy array.
            #Samples (column) by Taxa (rows)(lists) without the column
            data = xData.funcToArray()
            if data==None:
                print("Ordination:loadData::Error when converting AbundanceTable to Array, did not perform ordination.")
                return False

            #Transpose data to be Taxa (columns) by samples (rows)(lists)
            data = UtilityMath.funcTransposeDataMatrix(data,fRemoveAdornments=False)
            if(ValidateData.funcIsFalse(data)):
                print("Ordination:loadData::Error when transposing data file, did not perform ordination.")
                return False
            else:
                self.dataMatrix=data
                self.isRawData=fIsRawData
                self.lsIDs=xData.funcGetMetadata(xData.funcGetIDMetadataName())

        #Otherwise load the data directly as passed.
        else:
            self.dataMatrix=xData
            self.isRawData=fIsRawData
        return True

