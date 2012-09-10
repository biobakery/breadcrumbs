"""
Author: Timothy Tickle
Description: Plots matrices.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
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
