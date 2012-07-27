#######################################################
#
#	Title:		Constants_MicropitaTest
#	Author:		Timothy Tickle
#	Date:		September 18, 2011
#	Purpose:	Constants associated with automated unit testing
#
#######################################################

##
#Holds global configuration constants
class Constants_FiguresBreadCrumbs():

    #Figure oriented
    c_strBackgroundColorName = "Invisible"
    c_strBackgroundColor = "255,255,255"
    c_strBackgroundColorWord = "white"
    c_strBackgroundColorLetter = "w"

    #PCOA Markers
    c_charPCOAPieChart = "o"
    c_charPCOASquarePieChart = "s"
    iMarkerSize = 100

    c_fInverted = False

    def invertColors(self,fInvert):
        if fInvert==True:
            #General colors
            self.c_strBackgroundColor = "0,0,0"
            self.c_strBackgroundColorTuple = (0,0,0)
            self.c_strBackgroundColorWord = "black"
            self.c_strBackgroundColorLetter = "k"

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
            self.c_strDetailsColorWord = "white"
            self.c_strDetailsColorLetter = "w"

            #No select color
            self.c_charNoSelect = "#FFFFFF" # White

            #Record that it is not inverted
            self.c_fInverted = False

            #Alpha looks best at full in inversion
            self.c_dAlpha = 0.5
