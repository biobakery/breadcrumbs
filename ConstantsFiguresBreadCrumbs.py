"""
Author: Timothy Tickle
Description: Constants for figures.
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
#Holds global configuration constants
class ConstantsFiguresBreadCrumbs():

    #Figure oriented
    c_strBackgroundColorName = "Invisible"
    c_strBackgroundColor = "255,255,255"
    c_strBackgroundColorWord = "white"
    c_strBackgroundColorLetter = "w"
    c_strDetailsColorWord = "black"
    c_strDetailsColorLetter = "k"

    #PCOA Markers
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
