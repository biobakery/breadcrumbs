"""
Author: Timothy Tickle
Description: Class to create scatter plots.
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
class Histogram:

  @staticmethod
  def funcPlot(lx, strOutputFigurePath, strTitle = "Title", strXTitle="X Axis", strYTitle="Y Axis", strColor = "#83C8F9", fInvert=False):
    """
    Plot a box plot with optional jittering.

    :params	lx: List of x values
    :type:	List of doubles
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
    :params	fInvert: Invert colors (true)
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

    #Make scatter plot
    plt.hist(x=lx,histtype='bar',color=strColor)
    #Set ticks and title
    imgSubplot.set_title(strTitle)
    imgSubplot.title.set_color(objFigureControl.c_strDetailsColorLetter)

    #End plot
    #Save to a file
    imgFigure.savefig(strOutputFigurePath, facecolor=imgFigure.get_facecolor())
