"""
Author: Timothy Tickle
Description: Allows KMedoids on a custom metric space.
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

#External libraries
from scipy.spatial.distance import squareform

class MLPYDistanceAdaptor:
    """
    Allows one to use custom distance metrics with KMedoids in the MLPY package.
    """

    npaMatrix = None
    """
    Distance matrix to reference.
    """

    def __init__(self, npaDistanceMatrix, fIsCondensedMatrix):
        """
        Constructor requires a matrix of distances, could be condensed or square matrices

    	:param	npaDistanceMatrix:	The distance matrix to be used
	:type	Numpy array
	:param	fIsCondensedMatrix:	Indicator of the matrix being square (true = condensed; false = square)
	:type	Boolean
        """

        if fIsCondensedMatrix:
            self.npaMatrix = squareform(npaDistanceMatrix)
        else:
            self.npaMatrix = npaDistanceMatrix

    def compute(self,x,y):
        """
        This is the only method required in the interface to MLPY to be a distance metric.
        Does NOT want values but positions, the positions will be used for accessing the distance matrix already provided.

	:param	x:	X position as a array of 1 number
	:type	Numpy array
	:param	y:	Y position as a array of 1 number
	:type	Boolean
        """

        if(self.npaMatrix == None):
            raise Exception("".join(["MLPYDistanceAdaptor. Attempted to compute distance with out a distance matrix passed in during construction."]))
        return self.npaMatrix[x[0],y[0]]
