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
import os
from random import shuffle
from ValidateData import ValidateData

class SVM:
    """
    Class which holds generic methods for SVM use.
    """

    #1 Happy Path tested
    @staticmethod
    def funcConvertAbundanceTableToSVMFile(abndAbundanceTable, strOutputSVMFile, sMetadataLabel):
        """
        Converts abundance files to input SVM files.

        :param abndAbudnanceTable    AbudanceTable object to turn to input SVM file.
        :type    AbundnanceTable
        :param strOutputSVMFile: File to save SVM data to when converted from the abundance table.
        :type	String
        :param	sMetadataLabel: The name of the last row in the abundance table representing metadata.
        :type	String
        :return	lsUniqueLabels:	List of unique labels.
        """

        #Validate parameters
        if abndAbundanceTable == None:
            print "SVM::Error, invalid Abundance table."
            return False
        if(not ValidateData.funcIsValidString(strOutputSVMFile)):
            print "SVM::Error, file not valid. File:"+str(strOutputSVMFile)
            return False

        #If output file exists, delete
        if(os.path.exists(strOutputSVMFile)):
            os.remove(strOutputSVMFile)

        #Create data matrix
        dataMatrix = zip(*abndAbundanceTable.funcGetAbundanceCopy())

        #Add labels
        llData = []
        lsLabels = abndAbundanceTable.funcGetMetadata(sMetadataLabel)

        #Do not use set to make elements unique. Need to preserve order.
        #First label should be 0
        lsUniqueLabels = []
        for sElement in lsLabels:
            if sElement not in lsUniqueLabels:
                lsUniqueLabels.append(sElement)

        dictLabels = dict([[str(lenuLabels[1]),str(lenuLabels[0])] for lenuLabels in enumerate(lsUniqueLabels)])
        lsLabels = [dictLabels[sLabel] for sLabel in lsLabels]

        iRowIndex = 0
        for dataRow in dataMatrix[1:]:
            llData.append(ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace.join([lsLabels[iRowIndex]]+[ConstantsBreadCrumbs.c_strColon.join([str(enuSamples[0]+1),str(enuSamples[1])])
                            for enuSamples in enumerate(dataRow)])+ConstantsBreadCrumbs.c_strEndline)
            iRowIndex = iRowIndex + 1

        #Output file
        with open(strOutputSVMFile,'a') as f:
            (f.write("".join(llData)))

        return lsUniqueLabels

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
        for sElement in lsMetadata:
            if sElement not in lsUniqueLabels:
                lsUniqueLabels.append(sElement)

        dictLabels = dict([[str(lenuLabels[1]),str(lenuLabels[0])] for lenuLabels in enumerate(lsUniqueLabels)])
        return [dictLabels[sLabel] for sLabel in lsMetadata]

    #Tested
    @staticmethod
    def funcReadLabelsFromFile(sSVMFile, lsAllSampleNames, isPredictFile):
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

      :param sSVMFile:  File path to read in prediction labels.
      :type String
      :param lsAllSampleNames List of sample ids in the order of the labels.
      :type List of Strings
      :param isPredictFile: Indicates if the file is the input (False) or prediction (True) file
      :type boolean
      :return: Dictionary {label:["sampleName1", "sampleName2"...],...} or False on error
      """
      dictSampleLabels = dict()
 
      #Open prediction file and input file and get labels to compare to the predictions
      with open(sSVMFile,'r') as g:
          lsFileData = g.read()
          lsFileRows = [strRow for strRow in lsFileData.split(ConstantsBreadCrumbs.c_strEndline) if not(strRow is None) and not(strRow.strip() == '')]
          lsOriginalLabels = [row.split(ConstantsBreadCrumbs.c_strBreadCrumbsSVMSpace)[0] for row in lsFileRows[0+isPredictFile:]]
          #Check sample name length
          if not len(lsAllSampleNames) == len(lsOriginalLabels):
              print "SVM::Error, the length of sample names did not match the original labels length."
              return False
          for sValue in set(lsOriginalLabels):  
              dictSampleLabels[sValue] = set([lsAllSampleNames[iindex] for iindex, sLabel in enumerate(lsOriginalLabels) if sLabel == sValue])
      return dictSampleLabels

    #Tested
    @staticmethod
    def funcScaleFeature(npdData):
        """
        Scale a feature between 0 and 1. Using 01 and not 01,1 because it keeps te sparsity of the data and may save time.

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
