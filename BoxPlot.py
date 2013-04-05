"""
Author: Timothy Tickle
Description: Class to create box plots.
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
from ConstantsFiguresBreadCrumbs import ConstantsFiguresBreadCrumbs
import matplotlib.pyplot as plt
from pylab import *

#Plots a matrix
class BoxPlot:

  @staticmethod
  def funcPlot(ly, lsLabels, strOutputFigurePath, strTitle = "Title", strXTitle="X Axis", strYTitle="Y Axis", strColor = "#83C8F9", fJitter=False, fInvert=False, fInvertY=False):
    """
    Plot a box plot with optional jittering.

    :params	ly: List of y values
    :type:	List of doubles
    :params	lsLabels: List of labels (x tick lables)
    :type:	List string
    :params	strOutputFigurePath: File path to make figure
    :type:	String file path
    :params	strTitle: Title of figure
    :type:	String
    :params	strXTitle: Label of x axis
    :type:	String
    :params	strYTitle: Label of y axis
    :type:	String
    :params	strColor: Hex color for the face of the boxplots
    :type:	String
    :params	fJitter: Indicator of jittering (true) or not (false)
    :type:	Boolean
    :params	fInvert: Invert colors (true)
    :type:	Boolean
    :params	fInvertY: Invert y axis
    :type:	Boolean
    """

    #Start plot
    #Get plot object
    imgFigure = plt.figure()

    #Get plot colorsstrOutFigure
    objFigureControl = ConstantsFiguresBreadCrumbs()
    #Boxplots have to be plotted over the scatter so the alpha can not go to 1.0
    #In this case capturing the alpha before inversion
    #Inversion automoatically sets it to 1.
    dAlpha=objFigureControl.c_dAlpha
    objFigureControl.invertColors(fInvert=fInvert)

    #Color/Invert figure
    imgFigure.set_facecolor(objFigureControl.c_strBackgroundColorWord)
    imgSubplot = imgFigure.add_subplot(111,axisbg=objFigureControl.c_strBackgroundColorLetter)
    imgSubplot.set_xlabel(strXTitle)
    imgSubplot.set_ylabel(strYTitle)
    imgSubplot.spines['top'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.spines['bottom'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.spines['left'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.spines['right'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.xaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)

    #Adds light grid for numbers and puts them in the background
    imgSubplot.yaxis.grid(True, linestyle='-', which='major', color=objFigureControl.c_strGridLineColor, alpha=objFigureControl.c_dAlpha)
    imgSubplot.set_axisbelow(True)
    imgSubplot.yaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.tick_params(axis='x', colors=objFigureControl.c_strDetailsColorLetter)
    imgSubplot.tick_params(axis='y', colors=objFigureControl.c_strDetailsColorLetter)
    charMarkerEdgeColor = objFigureControl.c_strDetailsColorLetter

    #Make box plot
    bp = plt.boxplot(x=ly, notch=1, patch_artist=True)
    for iindex, ldData in enumerate(ly):
      ldX = None
      if fJitter:
        ldX = [float(iindex+1)+ uniform(-.05,.05) for x in xrange(len(ldData))]
      else:
        ldX = [float(iindex+1) for x in xrange(len(ldData))]
      plt.scatter(x=ldX,y=ly[iindex],c=strColor,marker="o",alpha=objFigureControl.c_dAlpha)

    #Color boxes
    plt.setp(bp['boxes'], color=objFigureControl.c_strDetailsColorLetter, facecolor=strColor, alpha=dAlpha)
    plt.setp(bp['whiskers'], color=objFigureControl.c_strDetailsColorLetter)

    #Set ticks and title
    lsLabelsWithCounts = []
    for iindex,sCurLabel in enumerate(lsLabels):
      lsLabelsWithCounts.append(sCurLabel+" ( "+str(len(ly[iindex]))+" )")
    xtickNames = plt.setp(imgSubplot, xticklabels=lsLabelsWithCounts)
    imgSubplot.set_title(strTitle)
    imgSubplot.title.set_color(objFigureControl.c_strDetailsColorLetter)

    #Invert Y axis
    if fInvertY:
      ax = plt.gca()
      ax.set_ylim(ax.get_ylim()[::-1])

    #End plot
    #Save to a file
    imgFigure.savefig(strOutputFigurePath, facecolor=imgFigure.get_facecolor())
