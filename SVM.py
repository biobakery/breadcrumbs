"""
Author: Timothy Tickle
Description: Class to Allow Support Vector Machine analysis and to contain associated scripts
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Libraries
from AbundanceTable import AbundanceTable
from ConstantsBreadCrumbs import ConstantsBreadCrumbs
import csv
import os
from random import shuffle
from ValidateData import ValidateData

class SVM:
    """
    Class which holds generic methods for SVM use.
    """

    #1 Happy Path tested
    @staticmethod
    def funcConvertAbundanceTableToSVMFile(abndAbundanceTable, xOutputSVMFile, sMetadataLabel, lsSampleOrdering = None):
        """
        Converts abundance files to input SVM files.

        :param abndAbundanceTable    AbudanceTable object to turn to input SVM file.
        :type    AbundanceTable
        :param strOutputSVMFile: File to save SVM data to when converted from the abundance table.
        :type	String
        :param	sMetadataLabel: The name of the last row in the abundance table representing metadata.
        :type	String
        :return	lsUniqueLabels:	List of unique labels.
        """

        #Create data matrix
        dataMatrix = zip(*abndAbundanceTable.funcGetAbundanceCopy())

        #Add labels
        llData = []
        lsLabels = SVM.funcMakeLabels(abndAbundanceTable.funcGetMetadata(sMetadataLabel))
        if not isinstance(xOutputSVMFile,str):
            if xOutputSVMFile.closed:
                xOutputSVMFile = open(xOutputSVMFile.name,"w")
	ostm = open(xOutputSVMFile,"w") if isinstance(xOutputSVMFile, str) else xOutputSVMFile
        f = csv.writer(ostm, csv.excel_tab, delimiter = ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)

	#This allows the creation of partially known files for stratification purposes
	lsCurrentSamples = abndAbundanceTable.funcGetSampleNames()
        lsOrderingSamples = lsSampleOrdering if lsSampleOrdering else lsCurrentSamples[:]

	iLabelIndex = 0
	iSize = len(dataMatrix[0])
	iIndexSample = 1
	for sSample in lsOrderingSamples:
		if sSample in lsCurrentSamples:
        		f.writerow([lsLabels[iLabelIndex]]+
				[ConstantsBreadCrumbs.c_strColon.join([str(tpleFeature[0]+1),str(tpleFeature[1])]) for tpleFeature in enumerate(dataMatrix[iIndexSample])])
			iLabelIndex += 1
			iIndexSample += 1
		#Make blank entry
		else:
			f.writerow([ConstantsBreadCrumbs.c_strSVMNoSample]+[ConstantsBreadCrumbs.c_strColon.join([str(tpleNas[0]+1),str(tpleNas[1])])
						for tpleNas in enumerate([ConstantsBreadCrumbs.c_strSVMNoSample]*iSize)])
	ostm.close()
        return set(lsLabels)

    @staticmethod
    def funcUpdateSVMFileWithAbundanceTable(abndAbundanceTable, xOutputSVMFile, lsOriginalLabels, lsSampleOrdering):
        """
        Converts abundance files to input SVM files.

        :param abndAbundanceTable    AbudanceTable object to turn to input SVM file.
        :type    AbundanceTable
        :param strOutputSVMFile: File to save SVM data to when converted from the abundance table.
        :type	String
        :param	sMetadataLabel: The name of the last row in the abundance table representing metadata.
        :type	String
        :return	lsUniqueLabels:	List of unique labels.
        """

        #Read in old file
        if not isinstance(xOutputSVMFile,str):
            if xOutputSVMFile.closed:
                xOutputSVMFile = open(xOutputSVMFile.name,"r")
	ostm = open(xOutputSVMFile,"r") if isinstance(xOutputSVMFile, str) else xOutputSVMFile
        fin = csv.reader(ostm, csv.excel_tab, delimiter = ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)
	#Read in contents of file
	llsOldContents = [lsRow for lsRow in fin]
	ostm.close()

	#Check to make sure this ordering covers all positions in the old file
	if not len(llsOldContents) == len(lsSampleOrdering):
		print "The length of the original file ("+str(len(llsOldContents))+") does not match the length of the ordering given ("+str(len(lsSampleOrdering))+")."
		return False

        #Create data matrix from new data
        dataMatrix = zip(*abndAbundanceTable.funcGetAbundanceCopy())

        #Add labels
        llData = []

	#Write to file
        if not isinstance(xOutputSVMFile,str):
            if xOutputSVMFile.closed:
                xOutputSVMFile = open(xOutputSVMFile.name,"w")
	ostm = open(xOutputSVMFile,"w") if isinstance(xOutputSVMFile, str) else xOutputSVMFile
        f = csv.writer(ostm, csv.excel_tab, delimiter = ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)

	#This allows to know what position to place te new lines
	lsCurrentSamples = abndAbundanceTable.funcGetSampleNames()

	iSize = len(dataMatrix[0])
	iIndexSample = 1
	iIndexOriginalOrder = 0
	for sSample in lsSampleOrdering:
		if sSample in lsCurrentSamples:
        		f.writerow([lsOriginalLabels[iIndexOriginalOrder]]+
				[ConstantsBreadCrumbs.c_strColon.join([str(tpleFeature[0]+1),str(tpleFeature[1])]) for tpleFeature in enumerate(dataMatrix[iIndexSample])])
			iIndexSample += 1
		#Make blank entry
		else:
			f.writerow(llsOldContents[iIndexOriginalOrder])
		iIndexOriginalOrder += 1
	ostm.close()
        return True

    #Tested 5
    @staticmethod
    def funcMakeLabels(lsMetadata):
        """
        Given a list of metadata, labels are assigned. This is function represents a central location to make labels so all are consistent.

        :param lsMetafdata:    List of metadata to turn into labels based on the metadata's values.
        :type    List of labels:    List of integer labels
        """
        #Do not use a set to make elements unique. Need to preserve order.
        #First label should be 0
        lsUniqueLabels = []
        [lsUniqueLabels.append(sElement) for sElement in lsMetadata if not (sElement in lsUniqueLabels)]

        dictLabels = dict([[str(lenuLabels[1]),str(lenuLabels[0])] for lenuLabels in enumerate(lsUniqueLabels)])
        return [dictLabels[sLabel] for sLabel in lsMetadata]

    #Tested
    @staticmethod
    def funcReadLabelsFromFile(xSVMFile, lsAllSampleNames, isPredictFile):
      """
      Reads in the labels from the input file or prediction output file of a LibSVM formatted file
      and associates them in order with the given sample names.

      Prediction file expected format: Labels declared in first line with labels keyword.
      Each following row a sample with the first entry the predicted label
      Prediction file example:
      labels 0 1
      0	0.3	0.4	0.6
      1	0.1	0.2	0.3
      1	0.2	0.2	0.2
      0	0.2	0.4	0.3

      Input file expected format:
      Each row a sample with the first entry the predicted label
      Input file example:
      0	0.3	0.4	0.6
      1	0.1	0.2	0.3
      1	0.2	0.2	0.2
      0	0.2	0.4	0.3

      :param xSVMFile:  File path to read in prediction labels.
      :type String
      :param lsAllSampleNames List of sample ids in the order of the labels.
      :type List of Strings
      :param isPredictFile: Indicates if the file is the input (False) or prediction (True) file
      :type boolean
      :return: Dictionary {label:["sampleName1", "sampleName2"...],...} or False on error
      """
      #Open prediction file and input file and get labels to compare to the predictions
      g = csv.reader( open(xSVMFile, 'r') if isinstance(xSVMFile, str) else xSVMFile, csv.excel_tab, delimiter = ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace )
      lsOriginalLabels = [lsLineElements[0] for lsLineElements in g if not lsLineElements[0] == ConstantsBreadCrumbs.c_strSVMNoSample]

      if isPredictFile:
          lsOriginalLabels = lsOriginalLabels[1:]

      #Check sample name length
      if not len(lsAllSampleNames) == len(lsOriginalLabels):
        print "SVM::funcReadLabelsFromFile. Error, the length of sample names did not match the original labels length. Samples ("+str(len(lsAllSampleNames))+"):"+str(lsAllSampleNames)+" Labels ("+str(len(lsOriginalLabels))+"):"+str(lsOriginalLabels)
        return False

      #Change to {label:["sampleName1", "sampleName2"...],...}
      dictSampleLabelsRet = dict()
      for sValue in set(lsOriginalLabels):  
        dictSampleLabelsRet[sValue] = set([lsAllSampleNames[iindex] for iindex, sLabel in enumerate(lsOriginalLabels) if sLabel == sValue])
      return dictSampleLabelsRet

    #Tested
    @staticmethod
    def funcScaleFeature(npdData):
        """
        Scale a feature between 0 and 1. Using 01 and not 01,1 because it keeps the sparsity of the data and may save time.

        :param	npdData:	Feature data to scale.
        :type	Numpy Array	Scaled feature data.
        :return npaFloat:    A numpy array of floats.
        """
        if sum(npdData) == 0 or len(set(npdData))==1:
            return npdData
        dMin = min(npdData)
        return (npdData-dMin)/float(max(npdData-dMin))

    #Tested
    @staticmethod
    def funcWeightLabels(lLabels):
        """
        Returns weights for labels based on how balanced the labels are. Weights try to balance unbalanced results.

        :params	lLabels:	List of labels to use for measure how balanced the comparison is.
        :type	List
        :return	List:		[dictWeights ({"label":weight}),lUniqueLabels (unique occurences of original labels)]
        """
        #Convert to dict
        #Do not use set to make elements unique. Need to preserve order.
        #First label should be 0
        lUniqueLabels = []
        for sElement in lLabels:
            if sElement not in lUniqueLabels:
                lUniqueLabels.append(sElement)
        dictLabels = dict(zip(lUniqueLabels, range(len(lUniqueLabels))))

        #Build a dict of weights per label {label:weight, label:weight}
        #Get the occurrence of each label
        dictWeights = dict()
        for sLabelKey in dictLabels:
            sCurLabel = dictLabels[sLabelKey]
            dictWeights[sCurLabel] = lLabels.count(sLabelKey)

        #Divide the highest occurrence each occurrence
        iMaxOccurence = max(dictWeights.values())
        for sWeightKey in dictWeights:
            dictWeights[sWeightKey]=iMaxOccurence/float(dictWeights[sWeightKey])

        return [dictWeights,lUniqueLabels]

    #Tested 3/4 cases could add in test 12 with randomize True
    def func10FoldCrossvalidation(self, iTotalSampleCount, fRandomise = False):
        """
        Generator.
        Generates the indexes for a 10 fold cross validation given a sample count.
        If there are less than 10 samples, it uses the sample count as the K-fold cross validation
        as a leave one out method.

        :param	iTotalSampleCount:	Total Sample Count
	    :type	Integer	Sample Count
	    :param	fRandomise:	Random sample indices
	    :type	Boolean	True indicates randomise (Default False)
        """
        #Make indices and shuffle if needed
        liindices = range(iTotalSampleCount)
        if fRandomise:
            shuffle(liindices)

        #For 10 times
        iKFold = 10
        if iTotalSampleCount < iKFold:
            iKFold = iTotalSampleCount
        for iiteration in xrange(iKFold):
            lfTraining = [iindex % iKFold != iiteration for iindex in liindices]
            lfValidation = [not iindex for iindex in lfTraining]
            yield lfTraining, lfValidation
