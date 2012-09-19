"""
Author: Timothy Tickle
Description: Plots matrices.
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
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#Plots a matrix
class PlotMatrix:

  #Given a matrix and labels consistent to the matrix, plot a matrix
  @staticmethod
  def funcPlotMatrix(npMatrix, lsXLabels, strOutputFigurePath, strXTitle="X Axis", strYTitle="Y Axis", fFlipYLabels=False):
    """
    Given a matrix and labels consistent to the matrix, plot a matrix.

    :param npMatrix: Numpy Array (matrix) to plot.
    :type: Numpy Array
    :param lsXLabels: X Labels
    :type: List of strings
    :param strOutputFigurePath: File to create the figure file.
    :type: String
    :param strXTitle: X Axis label.
    :type: String
    :param strYTitle: Y axis label.
    :type: String
    :param fFlipYLabels: Flip the Y labels so they are opposite order of x axis.
    :type: Boolean
    """

    #Get canvas/figure
    plt.clf()
    figConfusionMatrix = plt.figure()
    objAxis = figConfusionMatrix.add_subplot(111)

    #Get y labels
    lNewYLabels = list(lsXLabels)
    if fFlipYLabels:
        lNewYLabels.reverse()

    #Set x axis and position
    objAxis.xaxis.set_ticklabels([""]+lsXLabels)
    objAxis.xaxis.set_ticks_position('top')

    #Set y axis
    objAxis.yaxis.set_ticklabels([""]+lNewYLabels)

    #Set axis titles
    ylabel(strYTitle)
    plt.suptitle(strXTitle)

    #Plot matrix values
    objPlot = objAxis.imshow(np.array(npMatrix), cmap=get_cmap("Blues"), interpolation='nearest')

    #Plot text values
    for yIndex, ldRow in enumerate(npMatrix):
        for xIndex, dValue in enumerate(ldRow):
            plt.text(xIndex, yIndex, dValue, fontdict = {'size':18,'weight':'bold'} )

    #Add color bar
    figConfusionMatrix.colorbar(objPlot, ticks=range(int(min(np.array(npMatrix).ravel())),int(max(np.array(npMatrix).ravel()))))

    #Save to a file
    savefig(strOutputFigurePath)
