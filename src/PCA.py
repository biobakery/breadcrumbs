"""
Author: Timothy Tickle
Description: Performs and plots Principle Components Analysis.
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
from AbundanceTable import AbundanceTable
from ConstantsFiguresBreadCrumbs import ConstantsFiguresBreadCrumbs
from Ordination import Ordination
import matplotlib.cm as cm
from math import sqrt,asin
from matplotlib.mlab import PCA as mplPCA
from matplotlib import pyplot as plt
from numpy import *
from UtilityMath import UtilityMath
from ValidateData import ValidateData

class PCA(Ordination):
  """
  Class to Run Principle Components Analysis on an abundance table object
  """

  def __init__(self):
    Ordination.__init__(self)
    self.c_strComponents = "components"
    self.c_strVariance = "percent_variance"

  def run(self,fScale=True,fCenter=True,fASTransform=False):
    if not self.dataMatrix is None:
      mtrxPrepped = self.dataMatrix.T
      if fASTransform:
        mtrxPrepped = array([self.doAsinOnList(row) for row in sqrt(mtrxPrepped)])
      if fCenter:
        mtrxPrepped = mtrxPrepped-mean(mtrxPrepped,0)
      if fScale:
        # This is consistent to R's prcomp method.
        vStd = std(a=mtrxPrepped,axis=0) if fCenter else [sqrt(sum(square(ldRow))/len(ldRow)) for ldRow in mtrxPrepped.T]
        mtrxPrepped /= vStd
      iRows, iCols = mtrxPrepped.shape
      U,S,V = linalg.svd(a=mtrxPrepped,full_matrices=False)
      ldVariance = square(S*(iCols-1))
      ldVariance = ldVariance/sum(ldVariance)
      # Here components are row-wise so each component is a row.
      # Here percent variance is given and it is in the order of the components.
      self.dataProcessed = {self.c_strComponents:V, self.c_strVariance:ldVariance}
      return True
    else:
      print("PCA:run::Error Tried to run analysis on no data load data first.")
    return False

  def getVariance(self,iIndex=None):
    if not self.dataProcessed is None:
      if not iIndex is None:
        return self.dataProcessed[self.c_strVariance][iIndex]
      return self.dataProcessed[self.c_strVariance]
    else:
      print("PCA:getVariance::Error Tried to run analysis on no data load data first.")
    return False

  def getComponents(self,iIndex=None):
    if not self.dataProcessed is None:
      if not iIndex is None:
        return self.dataProcessed[self.c_strComponents].T[iIndex]
      return self.dataProcessed[self.c_strComponents].T
    else:
      print("PCA:getComponents::Error Tried to run analysis on no data load data first.")
    return False

  def doAsinOnList(self, lsValues):
    return([asin(element) for element in lsValues])

  def convertMetadataForPCA(self,abndTable):
    """ This takes a metadata dictionary from an abundance table and formats the metadata for use in the PCA.
        This formatting includes reducing discontinuous data to leveles and replacing NA values to the means of the value (continuous data only)
        This returns a numpy array of the format needed for this PCA object.
    """

    # Replace missing values with the mean
    # dummy the discrete data
    dictMetadata = abndTable.funcGetMetadataCopy()
    if(len(dictMetadata) < 2):
      return None

    ## Remove the metadata id
    dictMetadata.pop(abndTable.funcGetIDMetadataName(),None)
    lMetadata = []
    for lxItem in dictMetadata.values():
      ## If this is not numeric data then dummy
      ## Treat NA as a seperate category
      if not (sum([ ValidateData.funcIsValidStringFloat(xItem) for xItem in lxItem]) == len(lxItem)):
        # Get levels
        setsLevels = set(lxItem)
        # Go through each level and dummy the metadata
        for sLevel in setsLevels:
          lMetadata.append([1.0 if xItem==sLevel else 0.0 for xItem in lxItem])
      else:
        # Change NA to Mean and store numeric data as float
        # Also add to the metadata so that there are no negative numbers
        ldNONA = [float(xItem) for xItem in lxItem if not xItem.strip().lower() in ["na",""]]
        dMean = sum(ldNONA)/float(len(ldNONA))
        lsMetadataValues = [dMean if xItem.strip().lower() in ["na",""] else float(xItem) for xItem in lxItem]
        dMinValueAdj = abs(min(lsMetadataValues))
        lMetadata.append([sValue + dMinValueAdj for sValue in lsMetadataValues])
    return(array(lMetadata).T)
