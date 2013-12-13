"""
Author: Timothy Tickle
Description: Constants for figures.
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
#Holds global configuration constants
class ConstantsFiguresBreadCrumbs():

    #Figure oriented
    c_strBackgroundColorName = "Invisible"
    c_strBackgroundColor = "255,255,255"
    c_strBackgroundColorWord = "white"
    c_strBackgroundColorLetter = "w"
    c_strDetailsColorWord = "black"
    c_strDetailsColorLetter = "k"

    #PCoA Markers
    c_charPCOAPieChart = "o"
    c_charPCOASquarePieChart = "s"
    iMarkerSize = 100

    #PCoA defaults
    c_strPCoALabelDefault = "Label"
    c_cPCoAColorDefault = 'g'
    c_cPCoAShapeDefault = 'o'
    c_cPCoASizeDefault = 20

    #General plotting
    c_strGridLineColor = "#CCCCCC"

    c_fInverted = False
    c_dAlpha = 0.5

    def invertColors(self,fInvert):
        if fInvert==True:
            #General colors
            self.c_strBackgroundColor = "0,0,0"
            self.c_strBackgroundColorTuple = (0,0,0)
            self.c_strBackgroundColorWord = "black"
            self.c_strBackgroundColorLetter = "k"
            self.c_strDetailsColorWord = "white"
            self.c_strDetailsColorLetter = "w"

            #Invert no select color
            self.c_charNoSelect = "#000000" # black

            #Record that it is inverted
            self.c_fInverted = True

            #Alpha looks best at full in inversion
            self.c_dAlpha = 1.0

        else:
            #General colors
            self.c_strBackgroundColor = "255,255,255"
            self.c_strBackgroundColorTuple = (255,255,255)
            self.c_strBackgroundColorWord = "white"
            self.c_strBackgroundColorLetter = "w"
            self.c_strDetailsColorWord = "black"
            self.c_strDetailsColorLetter = "k"

            #No select color
            self.c_charNoSelect = "#FFFFFF" # White

            #Record that it is not inverted
            self.c_fInverted = False

            #Alpha looks best at full in inversion
            self.c_dAlpha = 0.5
