"""
Author: Timothy Tickle
Description: Perfroms and plots Principle Coordinates Analysis.
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
from cogent.cluster.nmds import NMDS
import csv
import math
import matplotlib.cm as cm
from Metric import Metric
import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats.stats import spearmanr
from Utility import Utility
from UtilityMath import UtilityMath
from ValidateData import ValidateData
from matplotlib import pyplot as plt

class PCoA:
    """
    Class to Run Principle Coordinates Analysis.

    To run PCoA first load the AbundanceTable or distance matrix using the "load" method, 
    then use the "run" method to derive points, and then use "plot" to plot the graph.
    The process is structured in this way so that data is read once but can be transformed to different
    distance matricies and after analysis can be plotted with multiple sample highlighting.
    One can always reload or rerun data by calling the appropriate function.

    Supported beta diversity metrics include "braycurtis","canberra","chebyshev","cityblock","correlation",
	"cosine","euclidean","hamming","sqeuclidean",unifrac_unweighted","unifrac_weighted"
    """

    #Supported distance metrics
    c_BRAY_CURTIS="B_Curtis"
    c_SPEARMAN="spearman"

    #Holds the data Matrix
    dataMatrix=None
    #Indicates if the data matrix is raw data (True) or a distance matrix (False)
    isRawData=None
    # Holds current matrix ids
    lsIDs = None

    #Current pcoa object
    pcoa = None

    #Holds the most recently successful distance metric
    strRecentMetric = None

    #Current dimensions
    _iDimensions = 2

    #Get plot colors
    objFigureControl = ConstantsFiguresBreadCrumbs()

    #Forced X Axis
    ldForcedXAxis = None

    #Indices for the plot group dictionary
    c_iXPointIndex = 0
    c_iYPointIndex = 1
    c_iColorIndex = 2
    c_iMarkerIndex = 3
    c_iAlphaIndex = 4
    c_iLabelIndex = 5
    c_iShapeIndex = 6
    c_iEdgeColorIndex = 7
    c_strTiesKey = "Ties"

    #Happy path tested
    def loadData(self, xData, fIsRawData):
        """
        Loads data into PCoA (given the matrix or an abundance table)
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
                print("PCoA:loadData::Error when converting AbundanceTable to Array, did not perform PCoA.")
                return False

            #Transpose data to be Taxa (columns) by samples (rows)(lists)
            data = UtilityMath.funcTransposeDataMatrix(data,fRemoveAdornments=False)
            if(ValidateData.funcIsFalse(data)):
                print("PCoA:loadData::Error when transposing data file, did not perform PCoA.")
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

    def run(self, tempDistanceMetric=None, iDims=2, strDistanceMatrixFile=None, istrmTree=None, istrmEnvr=None):
        """
        Runs analysis on loaded data.

        :param tempDistanceMetric: The name of the distance metric to use when performing PCoA.
                                   None indicates a distance matrix was already given when loading and will be used.
                                   Supports "braycurtis","canberra","chebyshev","cityblock","correlation",
				   "cosine","euclidean","hamming","sqeuclidean",unifrac_unweighted","unifrac_weighted"
        :type: String Distance matrix name
        :param iDims: How many dimension to plot the PCoA graphs.
                      (This can be minimally 2; all combinations of dimensions are plotted).
                      iDims start with 1 (not index-based).
        :type: Integer Positive integer 2 or greater.
	:param strDistanceMatrixFile: If the underlying distance matrix should be output, this is the file to output to.
	:type: String Output file for distances of None for indicating it shoudl not be done.
	:param istrmTree: One of two files needed for unifrac calculations, this is the phylogeny of the features.
	:type: String Path to file
	:param istrmEnvr: One of two files needed for unifrac calculations, this is the environment file for the features.
	:type: String Path to file
        :return boolean: Indicator of success (True)
        """

        if iDims > 1:
            self._iDimensions = iDims

        #If distance metric is none, check to see if the matrix is a distance matrix
        #If so, run NMDS on the distance matrix
        #Otherwise return a false and do not run
        if(tempDistanceMetric==None):
            if(ValidateData.funcIsTrue(self.isRawData)):
                print("PCoA:run::Error, no distance metric was specified but the previous load was not of a distance matrix.")
                return False
            elif(ValidateData.funcIsFalse(self.isRawData)):
                self.pcoa = NMDS(dataMatrix, verbosity=0)
                return True
        
        #Make sure the distance metric was a valid string type
        if(not ValidateData.funcIsValidString(tempDistanceMetric)):
            print("PCoA:run::Error, distance metric was not a valid string type.")
            return False

        #Supported distances
	
        distanceMatrix = None
        if(tempDistanceMetric==self.c_SPEARMAN):
            distanceMatrix = Metric().funcGetDissimilarity(ldSampleTaxaAbundancies=self.dataMatrix, funcDistanceFunction=lambda u,v: spearmanr(u,v)[0])
        if(tempDistanceMetric in [Metric.c_strUnifracUnweighted,Metric.c_strUnifracWeighted]):
            distanceMatrix,lsLabels = Metric().funcGetBetaMetric(sMetric=tempDistanceMetric, istrmTree=istrmTree, istrmEnvr=istrmEnvr)
            self.lsIDs = lsLabels
        else:
            distanceMatrix = Metric().funcGetBetaMetric(npadAbundancies=self.dataMatrix, sMetric=tempDistanceMetric)
        if(ValidateData.funcIsFalse(distanceMatrix)):
            print "PCoA:run::Error, when generating distance matrix."
            return False

        # Make squareform
        distanceMatrix = squareform(distanceMatrix)

        # Writes distance measures if needed.
        if strDistanceMatrixFile:
            csvrDistance = csv.writer(open(strDistanceMatrixFile, 'w'))
            if self.lsIDs:
                csvrDistance.writerow(["ID"]+self.lsIDs)

            for x in xrange(distanceMatrix.shape[0]):
                strId = [self.lsIDs[x]] if self.lsIDs else []
                csvrDistance.writerow(strId+distanceMatrix[x].tolist())

        self.pcoa = NMDS(distanceMatrix, dimension=max(self._iDimensions,2), verbosity=0)
        self.strRecentMetric = tempDistanceMetric
        return True

    #TODO Test
    def funcGetCoordinates(self):
        return(self.pcoa.getPoints())

    #TODO Test
    def funcGetIDs(self):
        return(self.lsIDs)

    #Happy path tested
    def plot(self, tempPlotName="PCOA.png", tempColorGrouping=None, tempShape=None, tempLabels=None, tempShapeSize=None, tempAlpha = 1.0, tempLegendLocation="upper right", tempInvert=False, iDim1=1, iDim2=2, fPlotOutline=True):
        """
        Plots the provided data by the given distance matrix in the file.
        All lists should be in order in relation to each other.
 
        :param tempPlotName: Path of file to save figure.
        :type: String File path.
        :param tempColorGrouping: Colors for markers.
                                  If you want a marker with multiple colors (piewedges) for that marker give a list in the list of colors.
                                  For example ['r','r','r',['r','g','b']] This would make 3 red markers and 1 split into  3 wedges (red, green, and blue).
                                  This is only possible if you are using circle shapes ('o') or square shapes ('s').
        :type: Character or list of characters: Characters should be useable by matplotlib as a color.
        :param tempShape: Marker shapes. If you want to specify one shape for all markers then just pass a char/str for the marker not a list.
        :type: Character or list of characters. Characters should be useable by matplotlib as shapes.
        :param tempLabels: Labels associated with the coloring. Should be consistent with tempColorGrouping (both should be strings or lists of equal length).
        :type: String or list of Strings.
        :param tempShapeSize: Sizes of markers (points). If no list is given, all markers are given the same size.
        :type: Integer of list of integers:	1 or greater.
        :param tempAlpha: Value between 0.0 and 1.0 (0.0 being completely transparent, 1.0 being opaque).
        :type: Float 0.0-1.0.
        :param tempLegendLocation: Indicates where to put the legend.
        :type: String Either "upper right", "lower right", "upper left", "lower left".
        :param tempInvert: Allows the inverting of the figure.
        :type: boolean True inverts.
        :param iDim1: First dimension to plot.
        :type: Integer Greater than 1.
        :param iDim2: Second dimension to plot.
        :type: Integer Greater than 1.
        :param fPlotOutline: Draw outline line around markers
        :type: boolean (True indicates draw outline)
        :return boolean: Indicator of success (True)
        """

        if(not self.pcoa == None):

            #Get point count
            iDimensionOne = max(0,min(self._iDimensions-2, iDim1-1))
            iDimensionTwo = max(1,min(self._iDimensions-1, iDim2-1))
            adPoints = self.pcoa.getPoints()

            #This is 1-stress which is the amount of variance not explained by all dimensions
            #There is no precent variance, so I am trying this as a substitute
            dPercentVariance = int((1.0-self.pcoa.getStress())*100)
            ldXPoints = list(adPoints[:,iDimensionOne])
            if not (self.ldForcedXAxis == None):
                ldXPoints = self.ldForcedXAxis
            ldYPoints = list(adPoints[:,iDimensionTwo])
            iPointCount = len(ldXPoints)

            #Get plot object
            imgFigure = plt.figure()
            self.objFigureControl.invertColors(fInvert=tempInvert)

            #Manage Labels
            if tempLabels is None:
                tempLabels = [self.objFigureControl.c_strPCoALabelDefault] * iPointCount
            elif(ValidateData.funcIsValidList(tempLabels)):
              if not len(tempLabels) == iPointCount:
                print "PCoA::plot:Error, the list of labels was given but was not the same length as the points so nothing was plotted."
                print "PCoA::plot:tempLabels=", tempLabels
                print "PCoA::plot:Label list length=", len(tempLabels) 
                print "PCoA::plot:iPointCount=", iPointCount
                return False
            elif ValidateData.funcIsValidString(tempLabels):
                tempLabels = [tempLabels] * iPointCount
            else:
                print "PCoA::plot:tempLabels was of an unexpected type. Expecting None, List, string, or char."
                print tempLabels
                return False

            #Manage Colors
            if tempColorGrouping is None:
                tempColorGrouping = [self.objFigureControl.c_cPCoAColorDefault] * iPointCount
            elif(ValidateData.funcIsValidList(tempColorGrouping)):
              if not len(tempColorGrouping) == iPointCount:
                print "PCoA::plot:Error, the list of colors was given but was not the same length as the points so nothing was plotted."
                print "PCoA::plot:tempColorGrouping=", tempColorGrouping
                print "PCoA::plot:Color list length=", len(tempColorGrouping) 
                print "PCoA::plot:iPointCount=", iPointCount
                return False
            elif ValidateData.funcIsValidString(tempColorGrouping):
                tempColorGrouping = [tempColorGrouping] * iPointCount
            else:
                print "PCoA::plot:tempColorGrouping was of an unexpected type. Expecting None, List, string, or char."
                print tempColorGrouping
                return False

            #Manage tempShape
            if tempShape is None:
                tempShape = [self.objFigureControl.c_cPCoAShapeDefault] * iPointCount
            elif(ValidateData.funcIsValidList(tempShape)):
              if not len(tempShape) == iPointCount:
                print "PCoA::plot:Error, the list of shapes was given but was not the same length as the points so nothing was plotted."
                print "PCoA::plot:tempShape=", tempShape
                print "PCoA::plot:Shape list length=", len(tempShape) 
                print "PCoA::plot:iPointCount=", iPointCount
                return False
            elif ValidateData.funcIsValidString(tempShape):
                tempShape = [tempShape] * iPointCount
            else:
                print("PCoA::plot:tempShape was of an unexpected type. Expecting None, List, string, or char.")
                print tempShape
                return False

            #Manage tempShapeSize
            if tempShapeSize is None:
                tempShapeSize = [self.objFigureControl.c_cPCoASizeDefault] * iPointCount
            elif(ValidateData.funcIsValidList(tempShapeSize)):
              if not len(tempShapeSize) == iPointCount:
                print "PCoA::plot:Error, the list of sizes was given but was not the same length as the points so nothing was plotted."
                print "PCoA::plot:tempShapeSize=", tempShapeSize
                print "PCoA::plot:Size list length=", len(tempShapeSize) 
                print "PCoA::plot:iPointCount=", iPointCount
                return False
            elif(ValidateData.funcIsValidInteger(tempShapeSize)):
                tempShapeSize = [tempShapeSize] * iPointCount
            else:
                print "PCoA::plot:tempShapeSize was of an unexpected type. Expecting None, List, string, or char."
                print tempShapeSize
                return False

            #Color/Invert figure
            imgFigure.set_facecolor(self.objFigureControl.c_strBackgroundColorWord)
            imgSubplot = imgFigure.add_subplot(111,axisbg=self.objFigureControl.c_strBackgroundColorLetter)
            imgSubplot.set_xlabel("Dimension "+str(iDimensionOne+1)+" (1-Stress = "+str(dPercentVariance)+"% )")
            imgSubplot.set_ylabel("Dimension "+str(iDimensionTwo+1))
            imgSubplot.spines['top'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.spines['bottom'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.spines['left'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.spines['right'].set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.xaxis.label.set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.yaxis.label.set_color(self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.tick_params(axis='x', colors=self.objFigureControl.c_strDetailsColorLetter)
            imgSubplot.tick_params(axis='y', colors=self.objFigureControl.c_strDetailsColorLetter)
            charMarkerEdgeColor = self.objFigureControl.c_strDetailsColorLetter if fPlotOutline else "none"  

            #If given a list of colors, each color will be plotted individually stratified by shape
            #Plot colors seperately so the legend will pick up on the labels and make a legend
            if(ValidateData.funcIsValidList(tempColorGrouping)):
                if len(tempColorGrouping) == iPointCount:

                    #Dictionary to hold plotting groups
                    #Logistical to plot points as layers in an intelligent fashion
                    #{CountofPoints: [[plot info list]]} The list happends so ties can occur in the key
                    dictPlotGroups = dict()
 
                    #Check for lists in the list which indicate the need to plot pie charts
                    lfAreLists = [ValidateData.funcIsValidList(objColor) for objIndex, objColor in enumerate(tempColorGrouping)]

                    #Pie chart data seperated out
                    lsColorsPieCharts = None
                    lcShapesPieCharts = None
                    lsLabelsPieCharts = None
                    lsSizesPieCharts = None
                    ldXPointsPieCharts = None
                    ldYPointsPieCharts = None

                    #Split out piechart data
                    if sum(lfAreLists) > 0:
                        #Get lists of index that are and are not lists
                        liAreLists = []
                        liAreNotLists = []
                        curIndex = 0
                        for fIsList in lfAreLists:
                            if fIsList: liAreLists.append(curIndex)
                            else: liAreNotLists.append(curIndex)
                            curIndex = curIndex + 1

                        lsColorsPieCharts = Utility.reduceList(tempColorGrouping, liAreLists)
                        tempColorGrouping = Utility.reduceList(tempColorGrouping, liAreNotLists)

                        #Split out shapes
                        lcShapesPieCharts = Utility.reduceList(tempShape, liAreLists)
                        tempShape = Utility.reduceList(tempShape, liAreNotLists)

                        #Split out labels
                        lsLabelsPieCharts = Utility.reduceList(tempLabels, liAreLists)
                        tempLabels = Utility.reduceList(tempLabels, liAreNotLists)

                        #Split out sizes
                        lsSizesPieCharts = Utility.reduceList(tempShapeSize, liAreLists)
                        tempShapeSize = Utility.reduceList(tempShapeSize, liAreNotLists)

                        #Split out xpoints
                        ldXPointsPieCharts = Utility.reduceList(ldXPoints, liAreLists)
                        ldXPoints = Utility.reduceList(ldXPoints, liAreNotLists)

                        #Split out ypoints
                        ldYPointsPieCharts = Utility.reduceList(ldYPoints, liAreLists)
                        ldYPoints = Utility.reduceList(ldYPoints, liAreNotLists)

                    #Get unique colors and plot each individually
                    acharUniqueColors = list(set(tempColorGrouping))
                    for iColorIndex in xrange(0,len(acharUniqueColors)):
                        #Get the color
                        charColor = acharUniqueColors[iColorIndex]

                        #Get indices of colors
                        aiColorPointPositions = Utility.getIndices(tempColorGrouping,charColor)

                        #Reduce the labels by color
                        acharLabelsByColor = Utility.reduceList(tempLabels,aiColorPointPositions)

                        #Reduces sizes to indices if a list
                        reducedSizes = tempShapeSize
                        #Reduce sizes if a list
                        if(ValidateData.funcIsValidList(reducedSizes)):
                          reducedSizes = Utility.reduceList(reducedSizes,aiColorPointPositions)

                        #Reduce to the current color grouping
                        aiXPoints = Utility.reduceList(ldXPoints,aiColorPointPositions)
                        aiYPoints = Utility.reduceList(ldYPoints,aiColorPointPositions)

                        #There are 3 options for shapes which are checked in this order.
                        #1. 1 shape character is given which is used for all markers
                        #2. A list is given of marker characters or lists of decimals which will be used to make pie chart markers
                        #This is handled after the rest this block of code
                        #3. A list of char are given each indicating the marker for a sample
                        #If the shapes are not a list plot
                        #Otherwise plot per shape per color (can not plot list of shapes in matplotlib)
                        reducedShapes = tempShape
                        if(not ValidateData.funcIsValidList(reducedShapes)):
                          reducedShapes = reducedShapes[0]
                          dictPlotGroups.setdefault(len(aiXPoints), []).append([aiXPoints,aiYPoints,[charColor],reducedShapes,tempAlpha,tempLabels[tempColorGrouping.index(charColor)],reducedSizes,charMarkerEdgeColor])
                        #Shapes are supplied as a list so plot each shape
                        else:
                          #Reduce to shapes of the current colors
                          reducedShapes = Utility.reduceList(reducedShapes,aiColorPointPositions)
                          acharReducedShapesElements = list(set(reducedShapes))
                          #If there are multiple shapes, plot seperately because one is not allowed to plot them as a list
                          for aCharShapeElement in acharReducedShapesElements:
                            #Get indices
                            aiShapeIndices = Utility.getIndices(reducedShapes,aCharShapeElement)
                            #Reduce label by shapes
                            strShapeLabel = Utility.reduceList(acharLabelsByColor,aiShapeIndices)
                            #Reduce sizes by shapes
                            strShapeSizes = reducedSizes
                            if ValidateData.funcIsValidList(reducedSizes):
                              strShapeSizes = Utility.reduceList(reducedSizes,aiShapeIndices)
                            #Get points per shape
                            aiXPointsPerShape = Utility.reduceList(aiXPoints,aiShapeIndices)
                            aiYPointsPerShape = Utility.reduceList(aiYPoints,aiShapeIndices)
                            #Get sizes per shape
                            #Reduce sizes if a list
                            reducedSizesPerShape = reducedSizes
                            if(ValidateData.funcIsValidList(reducedSizes)):
                              reducedSizesPerShape = Utility.reduceList(reducedSizes,aiShapeIndices)
                            #Put plot data in dict of lists for later plotting
                            #Separate out the background printing
                            dictPlotGroups.setdefault(len(aiXPointsPerShape), []).append([aiXPointsPerShape,aiYPointsPerShape,[charColor],aCharShapeElement,tempAlpha,strShapeLabel[0],strShapeSizes,charMarkerEdgeColor])

                    #Plot each color starting with largest color amount to smallest color amount so small groups will not be covered up by larger groups
                    #Plot other colors in increasing order
                    for sPlotGroupKey in sorted(list(dictPlotGroups.keys()), reverse=True):
                        lslsCurPlotGroup = dictPlotGroups[sPlotGroupKey]
                        #Plot
                        for lsGroup in lslsCurPlotGroup:
                            imgSubplot.scatter(lsGroup[self.c_iXPointIndex],
                                           lsGroup[self.c_iYPointIndex],
                                           c = lsGroup[self.c_iColorIndex],
                                           marker = lsGroup[self.c_iMarkerIndex],
                                           alpha = lsGroup[self.c_iAlphaIndex],
                                           label = lsGroup[self.c_iLabelIndex],
                                           s = lsGroup[self.c_iShapeIndex],
                                           edgecolor = lsGroup[self.c_iEdgeColorIndex])
 
                    #Plot pie charts
                    if not lsColorsPieCharts is None:
                        self.plotWithPieMarkers(imgSubplot=imgSubplot, aiXPoints=ldXPointsPieCharts, aiYPoints=ldYPointsPieCharts, dSize=lsSizesPieCharts, llColors=lsColorsPieCharts, lsLabels=lsLabelsPieCharts, lcShapes=lcShapesPieCharts, edgeColor=charMarkerEdgeColor, dAlpha=tempAlpha)

            objLegend = imgSubplot.legend(loc=tempLegendLocation, scatterpoints=1, prop={'size':10})

            #Invert legend
            if(tempInvert):
              if objLegend:
                objLegend.legendPatch.set_fc(self.objFigureControl.c_strBackgroundColorWord)
                objLegend.legendPatch.set_ec(self.objFigureControl.c_strDetailsColorLetter)
                plt.setp(objLegend.get_texts(),color=self.objFigureControl.c_strDetailsColorLetter)

            #Make legend background transparent
            if objLegend:
              objLegendFrame = objLegend.get_frame()
              objLegendFrame.set_alpha(self.objFigureControl.c_dAlpha)

            imgFigure.savefig(tempPlotName, facecolor=imgFigure.get_facecolor())
            return True

    #Indirectly tested
    def plotWithPieMarkers(self, imgSubplot, aiXPoints, aiYPoints, dSize, llColors, lsLabels, lcShapes, edgeColor, dAlpha):
        """
        The all lists should be in the same order

        :param imgSubPlot: Image to plot to
        :type: Image
        :param aiXPoints: List of X axis points (one element per color list)
        :type: List of Floats
        :param aiYPoints: List of X axis points (one element per color list)
        :type: List of Floats
        :param dSize: double or List of doubles (one element per color list)
        :type: List of Floats
        :param llColors: List of Lists of colors, one list of colors is for 1 piechart/multiply highlighted feature
                         Example ["red","blue","green"] for a marker with 3 sections.
        :type: List of strings
        :param lsLabels: List of labels  (one element per color list).
        :type: List of Floats
        :param lcShapes: Indicates which shape of a pie chart to use, currently supported 'o' and 's'  (one element per color list).
        :type: List of characters
        :param edgeColor: One color entry for the edge of the piechart.
        :type: List of characters
        :param dAlpha: Value between 0.0 and 1.0 (0.0 being completely transparent, 1.0 being opaque).
        :type: Float 0.0-1.0.
        """

        #Zip up points to pairs
        xyPoints = zip(aiXPoints,aiYPoints)
        #For each pair of points
        for iIndex,dXY in enumerate(xyPoints):
            ldWedges = []
            #Get colors
            lcurColors = llColors[iIndex]
            #Get pie cut shape
            cPieChartType = lcShapes[iIndex]
            if cPieChartType == ConstantsFiguresBreadCrumbs().c_charPCOAPieChart:
                ldWedges = self.makePieWedges(len(lcurColors),20)
            elif cPieChartType == ConstantsFiguresBreadCrumbs().c_charPCOASquarePieChart:
                ldWedges = self.makeSquarePieWedges(len(lcurColors))
            for iWedgeIndex,dWedge in enumerate(ldWedges):
                imgSubplot.scatter(x=dXY[0], y=dXY[1], marker=(dWedge,0), s=dSize[iIndex], label=lsLabels[iIndex], facecolor=lcurColors[iWedgeIndex], edgecolor=edgeColor, alpha=dAlpha)

    #Indirectly tested
    def makePieWedges(self, iWedgeCount, iSplineResolution = 10):
        """
        Generate a list of tuple points which will draw a square broken up into pie cuts.

        :param iWedgeCount: The number of piecuts in the square.
        :type: Integer Number greater than 1.
        :param iSplineResolution: The amount of smoothing to the circle's outer edge, the higher the number the more smooth.
        :type: integer Greater than 1.
        :return list List of tuples. Each tuple is a point, formatted for direct plotting of the marker.
        """

        ldWedge = []
        dLastValue = 0.0

        #Create a list of equal percentages for all wedges
        #Do not include a last wedge it gets all the space from the 2nd to last wedge to the end
        #Which should still be equal to the others
        ldPercentages = [1.0/iWedgeCount]*(iWedgeCount-1)

        for dPercentage in ldPercentages:
            ldX = [0] + np.cos(np.linspace(2*math.pi*dLastValue,2*math.pi*(dLastValue+dPercentage),iSplineResolution)).tolist()
            ldY = [0] + np.sin(np.linspace(2*math.pi*dLastValue,2*math.pi*(dLastValue+dPercentage),iSplineResolution)).tolist()
            ldWedge.append(zip(ldX,ldY))
            dLastValue = dLastValue+dPercentage
        ldX = [0] + np.cos(np.linspace(2*math.pi*dLastValue,2*math.pi,iSplineResolution)).tolist()
        ldY = [0] + np.sin(np.linspace(2*math.pi*dLastValue,2*math.pi,iSplineResolution)).tolist()
        ldWedge.append(zip(ldX,ldY))
        return ldWedge

    #Indirectly tested
    def makeSquarePieWedges(self, iWedgeCount):
        """
        Generate a list of tuple points which will draw a square broken up into pie cuts.

        :param iWedgeCount: The number of piecuts in the square.
        :type: Integer Number greater than 1.
        :return list List of tuples. Each tuple is a point, formatted for direct plotting of the marker.
        """

        ldWedge = []
        dLastPercentageValue = 0.0
        dLastSquareValue = 0.0
        dCumulativePercentageValue = 0.0
        dRadius = None
        fXYSwitched = False
        fAfterCorner = False
        iSwitchCounts = 0
        iMagicNumber =(1.0/4)

        #Create a list of equal percentages for all wedges
        #Do not include a last wedge it gets all the space from the 2nd to last wedge to the end
        #Which should still be equal to the others
        ldPercentages = [1.0/iWedgeCount]*(iWedgeCount)

        for dPercentage in ldPercentages:
          ldCircleXs = np.cos([2*math.pi*dLastPercentageValue,2*math.pi*(dLastPercentageValue+dPercentage)])
          ldCircleYs = np.sin([2*math.pi*dLastPercentageValue,2*math.pi*(dLastPercentageValue+dPercentage)])

          if dRadius == None:
            dRadius = ldCircleXs[0]

          #Check to see if at corner
          fAtCorner = False
          iDistance = int((dLastPercentageValue+dPercentage+(iMagicNumber/2))/iMagicNumber
                  ) - int((dLastPercentageValue+(iMagicNumber/2))/iMagicNumber)
          if(iDistance > 0):
            fAtCorner = True
            if iDistance > 1:
              fXYSwitched = not fXYSwitched
              iSwitchCounts = iSwitchCounts + 1

          #Check to see if at a side center
          fAtSide = False
          if (int((dLastPercentageValue+dPercentage)/iMagicNumber) > int(dLastPercentageValue/iMagicNumber)):
            fAtSide = True

          #Handle corner xy switching
          if fAtCorner:
            fXYSwitched = not fXYSwitched
            iSwitchCounts = iSwitchCounts + 1
          #Make sure the xy switching occurs to vary the slope at the corner.
          if fXYSwitched:
              ldCircleXs,ldCircleYs = ldCircleYs,ldCircleXs

          dSquarePoint = dRadius * (ldCircleYs[1]/float(ldCircleXs[1]))
          dRadiusSq1 = dRadius
          dRadiusSq2 = dRadius
          dLastSquareValueSq = dLastSquareValue
          dSquarePointSq = dSquarePoint

          #If in quadrants 2,3 make sign changes
          if iSwitchCounts in [2,3]:
            if iSwitchCounts == 2:
              dRadiusSq1 = dRadiusSq1 *-1
            elif iSwitchCounts == 3:
              dRadiusSq1 = dRadiusSq1 * -1
              dRadiusSq2 = dRadiusSq2 * -1
            dLastSquareValueSq = dLastSquareValueSq * -1.0
            dSquarePointSq = dSquarePointSq * -1.0

          if fAtCorner:
            #Corner 1
            if iSwitchCounts==1:
              ldWedge.append(zip([0,dRadiusSq1,dRadiusSq1,dSquarePointSq,0],[0,dLastSquareValueSq,dRadiusSq2,dRadiusSq2,0]))
            #Corner 2
            elif iSwitchCounts==2:
              if iDistance > 1:
                ldWedge.append(zip([0,-dRadiusSq1,-dRadiusSq1,dRadiusSq1,dRadiusSq1,0],[0,-dLastSquareValueSq,dRadiusSq2,dRadiusSq2,dSquarePointSq,0]))
              else:
                ldWedge.append(zip([0,-dLastSquareValueSq,dRadiusSq1,dRadiusSq1,0],[0,dRadiusSq2,dRadiusSq2,dSquarePointSq,0]))
            #Corner 3
            elif iSwitchCounts==3:
              if iDistance > 1:
                ldWedge.append(zip([0,-dLastSquareValueSq,dRadiusSq1,dRadiusSq1,dSquarePointSq,0],[0,-dRadiusSq2,-dRadiusSq2,dRadiusSq2,dRadiusSq2,0]))
              else:
                ldWedge.append(zip([0,dRadiusSq1,dRadiusSq1,dSquarePointSq,0],[0,dLastSquareValueSq,dRadiusSq2,dRadiusSq2,0]))
            #Corner 4
            elif iSwitchCounts==4:
              if iDistance > 1:
                ldWedge.append(zip([0,-dRadiusSq1,-dRadiusSq1,dRadiusSq1,dRadiusSq1,0],[0,-dLastSquareValueSq,-dRadiusSq2,-dRadiusSq2,dSquarePointSq,0]))
              else:
                ldWedge.append(zip([0,(-1.0*dLastSquareValueSq),dRadiusSq1,dRadiusSq1,0],[0,(-1.0*dRadiusSq2),(-1.0*dRadiusSq2),dSquarePointSq,0]))

            fAfterCorner = True
          else:
            if iSwitchCounts%2:
              ldWedge.append(zip([0,dLastSquareValueSq,dSquarePointSq,0],[0,dRadiusSq2,dRadiusSq2,0]))
            else:
              ldWedge.append(zip([0,dRadiusSq1,dRadiusSq1,0],[0,dLastSquareValueSq,dSquarePointSq,0]))

          dLastSquareValue = dSquarePoint
          dCumulativePercentageValue = dCumulativePercentageValue + dLastSquareValue
          dLastPercentageValue = dLastPercentageValue+dPercentage

        return ldWedge

    #Happy Path Tested
    def plotList(self, lsLabelList, strOutputFileName, iSize=20, dAlpha=1.0, charForceColor=None, charForceShape=None, fInvert=False, iDim1=1, iDim2=2, fPlotOutline=True, sLegendLocation="upper right"):
        """
        Convenience method used to plot data in the PCoA given a label list (which is in order of the underlying data).
        This is for the scenario where you do not care that the color or shape of the data will be as long as it varies
        with the label.
        This method does allow forcing color or shape to 1 character so that they do not vary with the label but are one value.
        This is helpful when you have a large number of labels to plot given the shapes in the PCoA are limited but not the coloring.

        :param lsLabelList: List of string labels which are in order of the data in the PCoA object (as the data was loaded the PCoA object).
        :type: List of strings
        :param strOutputFileName: File path to save figure.
        :type: String
        :param iSize: Size of marker. Default 20.
        :type: Integer
        :param dAlpha: Alpha for the markers. (0.0 tranparent, 1.0 opaque)
        :type: Double between 0.0 and 1.0
        :param charForceColor: Color to force the points to. (Must be understandable by matplotlib as a color [ie. 'k','m','c','r','g','b','y','w'])
        :type: Character
        :param charForceShape: Shape to force the points to. (Must be understandable by matplotlib as a shape [ie. 'o','s','^','v','<','>','8','p','h']), False makes all shapes a circle.
        :type: Character or False
        :param fInvert: Allows one to invert the background and plot details from white to black (True == background is black).
        :type: Boolean
        :param iDim1: The first dimension to plot
        :type: Integer starting at 1
        :param iDim2: The second dimension to plot
        :type: Integer starting at 2
        :return boolean: Indicator of success (True)
        """

        #Get uniqueValues for labels
        acharUniqueValues = list(set(lsLabelList))
        iCountUniqueValues = len(acharUniqueValues)

        #Set colors
        atupldLabelColors = None

        #Set shapes
        alLabelShapes = None
        if charForceShape == None:
            #Get shapes
            acharShapes = PCoA.getShapes(iCountUniqueValues)
            if len(acharShapes) == 0:
                return False
            #Make label shapes
            alLabelShapes = [ acharShapes[acharUniqueValues.index(sMetadata)] for sMetadata in lsLabelList ]
        elif charForceShape == False:
            alLabelShapes = [ self.objFigureControl.c_cPCoAShapeDefault ] * len(lsLabelList)
        else:
            alLabelShapes = charForceShape

        #If the coloring is not forced, color so it is based on the labels
        if charForceColor == None:
            #Get colors based on labels
            atupldColors = [ 
                Utility.RGBToHex(
                    cm.jet( float(iUniqueValueIndex)/float(iCountUniqueValues) )
                    ) 
                for iUniqueValueIndex in xrange(0,iCountUniqueValues)
                ]
            #Make label coloring
            atupldLabelColors = [ atupldColors[acharUniqueValues.index(sMetadata) ] for sMetadata in lsLabelList ]
        elif type( charForceColor ) is dict:
            atupldLabelColors = [ charForceColor.get(sMetadata,self.objFigureControl.c_cPCoAColorDefault) for sMetadata in lsLabelList ]
        #If the coloring is forced, color so it is based on the charForcedColor list
        elif(ValidateData.funcIsValidList(charForceColor)):
            atupldLabelColors = charForceColor[0]
            if not len(lsLabelList) == len(atupldLabelColors):
                print "PCoA::plotList:Error, label and forced color lengths were not the same."
                print "Labels"
                print lsLabelList
                print len(lsLabelList)
                print "Forced Colors"
                print charForceColor[0]
                print len(charForceColor[0])
                return False
            lsLabelList = [ "".join([charForceColor[1][iLabelIndex], "_", lsLabelList[iLabelIndex]]) for iLabelIndex in xrange(0,len(charForceColor[1]))]
        #If the color is forced but the color does not vary, color all markers are the same.
        else:
            atupldLabelColors = charForceColor

        #Call plot
        self.plot(tempPlotName=strOutputFileName, tempColorGrouping=atupldLabelColors, tempShape=alLabelShapes, tempLabels=lsLabelList, tempShapeSize = iSize, tempAlpha=dAlpha, tempLegendLocation=sLegendLocation , tempInvert = fInvert, iDim1=iDim1, iDim2=iDim2, fPlotOutline=fPlotOutline)

    def funcForceXAxis(self, dList):
        """
        Force the X axis to the given list.

        :param dList: List of values to force the x axis of the plot (floats).
        :type: List of floats
        """

        self.ldForcedXAxis = dList

    def funcUnforceXAxis(self):
        """
        Return the X axis to the values derived from the loaded data.
        """

        self.ldForcedXAxis = None

    #Happy Path Tested
    @staticmethod
    def getShapes(intShapeCount):
        """
        Returns a list of characters which are valid shapes for markers.

        :param intShapeCount: The number of shapes to return.
        :type: Integer (min 1, max 9)
        :return: A list of characters to use as markers. [] is returned on error
        """

        lsPointShapes = ['o','s','^','v','<','>','8','p','h']
        if intShapeCount > len(lsPointShapes):
            print("".join(["Error, PCoA.getShapes. Do not have enough shapes to give. Received request for ",str(intShapeCount)," shapes. Max available shape count is ",str(len(lsPointShapes)),"."]))
            return []
        return lsPointShapes[0:intShapeCount]
