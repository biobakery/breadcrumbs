"""
Author: Timothy Tickle
Description: Calculates Metrics.
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

#Update path
from ConstantsBreadCrumbs import ConstantsBreadCrumbs
import csv
import numpy as np
from types import *
from ValidateData import ValidateData

#External libraries
from cogent.maths.unifrac.fast_unifrac import fast_unifrac_file
import cogent.maths.stats.alpha_diversity
import scipy.spatial.distance

class Metric:
    """
    Performs ecological measurements.
    """

    #Diversity metrics Alpha
    c_strSimpsonDiversity = "SimpsonD"
    c_strInvSimpsonDiversity = "InSimpsonD"
    c_strChao1Diversity = "Chao1"

    #Diversity metrics Beta
    c_strBrayCurtisDissimilarity = "B_Curtis"
    c_strUnifracUnweighted = "U_unifrac"
    c_strUnifracWeighted = "W_unifrac"

    #Additive inverses of beta metrics
    c_strInvBrayCurtisDissimilarity = "InB_Curtis"

    #Richness
    c_strShannonRichness = "ShannonR"
    c_strObservedCount = "Observed_Count"

    #Different count-based alpha diversity metrics
    setAlphaDiversities = set(["observed_species","margalef","menhinick",
	"dominance","reciprocal_simpson","shannon","equitability","berger_parker_d",
	"mcintosh_d","brillouin_d","kempton_taylor_q","strong","fisher_alpha",
	"mcintosh_e","heip_e","simpson_e","robbins","michaelis_menten_fit","chao1","ACE"])

    #Different count-based beta diversity metrics
    setBetaDiversities = set(["braycurtis","canberra","chebyshev","cityblock",
	"correlation","cosine","euclidean","hamming","sqeuclidean"])

    #Tested 4
    @staticmethod
    def funcGetSimpsonsDiversityIndex(ldSampleTaxaAbundancies=None):
        """
        Calculates the Simpsons diversity index as defined as sum(Pi*Pi).
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:	List of measurements to calculate metric on (a sample).
        :type:	List of doubles
        :return	Double:	Diversity metric
        """

        #Calculate metric
        return sum((ldSampleTaxaAbundancies)*(ldSampleTaxaAbundancies))

    #Tested 4
    @staticmethod
    def funcGetInverseSimpsonsDiversityIndex(ldSampleTaxaAbundancies=None):
        """
        Calculates Inverse Simpsons diversity index 1/sum(Pi*Pi).
        This is multiplicative inverse which reverses the order of the simpsons diversity index.
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:	List of measurements to calculate metric on (a sample).
        :type:	List of doubles
        :return	Double:	Diversity metric
        """

        simpsons = Metric.funcGetSimpsonsDiversityIndex(ldSampleTaxaAbundancies)
        #If simpsons is false return false, else return inverse
        if not simpsons:
            return False
        return 1.0/simpsons

    #Tested 4
    @staticmethod
    def funcGetShannonRichnessIndex(ldSampleTaxaAbundancies=None):
        """
        Calculates the Shannon richness index.
        Note***: Assumes that the abundance measurements are already normalized by the total population N.
        If not normalized, include N in the parameter tempTotalN and it will be.
	This is in base exp(1) like the default R Vegan package. Cogent is by defaul in bits (base=2)
	Both options are here for your use. See Metric.funcGetAlphaDiversity() to access cogent

        :param	ldSampleTaxaAbundancies:	List of measurements to calculate metric on (a sample).
        :type:	List of doubles
        :return	Double:	Richness metric
        """

        #Calculate metric
        ldSampleTaxaAbundancies = ldSampleTaxaAbundancies[np.where(ldSampleTaxaAbundancies != 0)]
        tempIntermediateNumber = sum(ldSampleTaxaAbundancies*(np.log(ldSampleTaxaAbundancies)))
        if(tempIntermediateNumber == 0.0):
            return 0.0
        return -1 * tempIntermediateNumber

    #Test 3
    @staticmethod
    def funcGetChao1DiversityIndex(ldSampleTaxaAbundancies=None, fCorrectForBias=False):
        """
        Calculates the Chao1 diversity index.
        Note***: Not normalized by abundance.

        :param	ldSampleTaxaAbundancies:	List of measurements to calculate metric on (a sample).
        :type:	List of doubles
        :param	fCorrectForBias:	Indicator to use bias correction.
        :type:	Boolean	False indicates uncorrected for bias (uncorrected = Chao 1984, corrected = Chao 1987, Eq. 2)
        :return	Double:	Diversity metric
        """
        #If not counts return false
        if [num for num in ldSampleTaxaAbundancies if((num<1) and (not num==0))]: return False

        #Observed = total number of species observed in all samples pooled
        totalObservedSpecies = len(ldSampleTaxaAbundancies)-len(ldSampleTaxaAbundancies[ldSampleTaxaAbundancies == 0])

        #Singles = number of species that occur in exactly 1 sample
        singlesObserved = len(ldSampleTaxaAbundancies[ldSampleTaxaAbundancies == 1.0])

        #Doubles = number of species that occue in exactly 2 samples
        doublesObserved = len(ldSampleTaxaAbundancies[ldSampleTaxaAbundancies == 2.0])

        #If singles or doubles = 0, return observations so that a divided by zero error does not occur
        if((singlesObserved == 0) or (doublesObserved == 0)):
            return totalObservedSpecies

        #Calculate metric
        if fCorrectForBias:
            return cogent.maths.stats.alpha_diversity.chao1_bias_corrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)
        else:
            return cogent.maths.stats.alpha_diversity.chao1_uncorrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)

    #Test 3
    @staticmethod
    def funcGetObservedCount(ldSampleAbundances, dThreshold = 0.0):
        """
        Count how many bugs / features have a value of greater than 0 or the threshold given.
        Expects a vector of abundances.
        ****Do not normalize data if using the threshold.

        :param	ldSampleAbundances:	List of measurements to calculate metric on (a sample).
        :type:	List of doubles
        :param	dThreshold:	The lowest number the measurement can be to be counted as an observation.
        :type:	Double
        :return	Count:	Number of features observed in a sample.
        """

        return sum([1 for observation in ldSampleAbundances if observation > dThreshold])

    #Test Cases 6
    @staticmethod
    def funcGetAlphaDiversity(liCounts,strMetric):
        """
        Passes counts to cogent for an alpha diversity metric.
	setAlphaDiversities are the names supported

        :param	liCount:	List of counts to calculate metric on (a sample).
        :type:	List of ints
        :return	Diversity:	Double diversity metric.
        """

        return getattr(cogent.maths.stats.alpha_diversity,strMetric)(liCounts)

    #Happy path tested 1
    @staticmethod
    def funcGetDissimilarity(ldSampleTaxaAbundancies, funcDistanceFunction):
        """
        Calculates the distance between samples given a function.

        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:
        :type:	List of doubles
        :param	funcDistanceFunction: Distance function used to calculate distances
        :type:	Function
        :return	Double:	Dissimilarity metric
        """

        #Calculate metric
        try:
            return scipy.spatial.distance.pdist(ldSampleTaxaAbundancies, funcDistanceFunction)
        except ValueError as error:
            print "".join(["Metric.funcGetDissimilarity. Error=",str(error)])
            return False

    #Test case 1
    @staticmethod
    def funcGetDissimilarityByName(ldSampleTaxaAbundancies, strMetric):
        """
        Calculates beta-diversity metrics between lists of abundances
	setBetaDiversities are the names supported

        :param	ldSampleTaxaAbundancies:
        :type:	List of doubles
        :param	strMetric: Name of the distance function used to calculate distances
        :type:	String
        :return	list double:	Dissimilarity metrics between each sample
        """

        return scipy.spatial.distance.pdist(ldSampleTaxaAbundancies,strMetric)

    #Test 3
    @staticmethod
    def funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies):
        """
        Calculates the BrayCurtis Beta dissimilarity index.
        d(u,v)=sum(abs(row1-row2))/sum(row1+row2).
        This is scale invariant.
        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:
        :type:	List of doubles
        :return	Double Matrix:	Dissimilarity metric
        """

        #Calculate metric
        try:
            return scipy.spatial.distance.pdist(X=ldSampleTaxaAbundancies, metric='braycurtis')
        except ValueError as error:
            print "".join(["Metric.getBrayCurtisDissimilarity. Error=",str(error)])
            return False

    #Test 3
    @staticmethod
    def funcGetInverseBrayCurtisDissimilarity(ldSampleTaxaAbundancies):
        """
        Calculates 1 - the BrayCurtis Beta dissimilarity index.
        d(u,v)=1-(sum(abs(row1-row2))/sum(row1+row2)).
        This is scale invariant and ranges between 0 and 1.
        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:	An np.array of samples (rows) x measurements (columns) in which distance is measured between rows
        :type:	List	List of doubles
        :return	Double Matrix:	1 - Bray-Curtis dissimilarity.	
        """

        bcValue = Metric.funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies = ldSampleTaxaAbundancies)
        if not type(bcValue) is BooleanType:
            return 1.0-bcValue
        return False

    #Test cases 8
    @staticmethod
    def funcGetUnifracDistance(istrmTree,istrmEnvr,lsSampleOrder=None,fWeighted=True):
	"""
	Gets a unifrac distance from files/filestreams.

        :param	istrmTree:	File path or stream which is a Newick format file
        :type:	String of file stream
        :param	istrmEnvr:	File path or stream which is a Newick format file
        :type:	String of file stream
	"""
	npaDist, lsSampleNames = fast_unifrac_file(open(istrmTree,"r") if isinstance(istrmTree, str) else istrmTree,
			open(istrmEnvr,"r") if isinstance(istrmEnvr, str) else istrmEnvr, weighted=fWeighted).get("distance_matrix",False)

        #Was trying to avoid preallocating a matrix but if you only need a subset of the samples then it
        #is simpler to preallocate so this is what I am doing but making a condensed matrix and not a full matrix
        
        #Dictionary to translate the current order of the samples to what is expected if given an input order
        if lsSampleOrder:
            #{NewOrder:OriginalOrder} way to convert from old to new sample location
            dictTranslate = dict([[lsSampleOrder.index(sSampleName),lsSampleNames.index(sSampleName)] for sSampleName in lsSampleNames if sSampleName in lsSampleOrder])

            #Check to make sure all samples requested were found
            if not len(dictTranslate.keys()) == len(lsSampleOrder):
                print "Metric.funcGetUnifracDistance. Error= The some or all sample names given (lsSampleOrder) were not contained in the matrix."
                return False

            #Length of data
            iLengthOfData = len(lsSampleOrder)

            #Preallocate matrix and shuffle
            mtrxData = np.zeros(shape=(iLengthOfData,iLengthOfData))
            for x in xrange(iLengthOfData):
                for y in xrange(iLengthOfData):
                    mtrxData[x,y] = npaDist[dictTranslate[x],dictTranslate[y]]
            npaDist = mtrxData

            lsSampleNames = lsSampleOrder

        #If no sample order is given, condense the matrix and return
        return (scipy.spatial.distance.squareform(npaDist),lsSampleNames)


    #Test 7
    @staticmethod
    def funcGetAlphaMetric(ldAbundancies, strMetric):
        """
        Get alpha abundance of the metric for the vector.
        Note: Shannon is measured with base 2 ("shannon") or base exp(1) (Metric.c_strShannonRichness) depending which method is called.

        :param	ldAbundancies:	List of values to compute metric (a sample).
        :type:	List	List of doubles.
        :param	strMetric:	The metric to measure.
        :type:	String	Metric name (Use from constants above).
        :return	Double:	Metric specified by strMetric derived from ldAbundancies.
        """

        if(strMetric == Metric.c_strShannonRichness):
            return Metric.funcGetShannonRichnessIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric == Metric.c_strSimpsonDiversity):
            return Metric.funcGetSimpsonsDiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric == Metric.c_strInvSimpsonDiversity):
            return Metric.funcGetInverseSimpsonsDiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric == Metric.c_strObservedCount):
            return Metric.funcGetObservedCount(ldSampleAbundances=ldAbundancies)
        #Chao1 Needs NOT Normalized Abundance (Counts)
        elif(strMetric == Metric.c_strChao1Diversity):
            return Metric.funcGetChao1DiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        elif(strMetric in Metric.setAlphaDiversities):
            return Metric.funcGetAlphaDiversity(liCounts=ldAbundancies, strMetric=strMetric)
        else:
            return False

    #Test 5
    @staticmethod
    def funcBuildAlphaMetricsMatrix(npaSampleAbundance = None, lsSampleNames = None, lsDiversityMetricAlpha = None):
        """
        Build a matrix of alpha diversity metrics for each sample
        Row = metric, column = sample

        :param	npaSampleAbundance:	Observations (Taxa (row) x sample (column))
        :type:	Numpy Array
        :param	lsSampleNames:	List of sample names of samples to measure (do not include the taxa id column name or other column names which should not be read).
        :type:	List of strings	Strings being samples to measure from the npaSampleAbundance.
        :param	lsDiversityMetricAlpha:	List of diversity metrics to use in measuring.
        :type:	List of strings	Strings being metrics to derived from the indicated samples.
        :return	List of List of doubles:	Each internal list is a list of (floats) indicating a specific metric measurement method measuring multiple samples
            [[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
        """

        if not ValidateData.funcIsValidList(lsDiversityMetricAlpha):
            lsDiversityMetricAlpha = [lsDiversityMetricAlpha]

        #Get amount of metrics
        metricsCount = len(lsDiversityMetricAlpha)

        #Create return
        returnMetricsMatrixRet = [[] for index in lsDiversityMetricAlpha]

        #For each sample get all metrics
        #Place in list of lists
        #[[metric1-sample1, metric1-sample2, metric1-sample3],[metric1-sample1, metric1-sample2, metric1-sample3]]
        for sample in lsSampleNames:
            sampleAbundance = npaSampleAbundance[sample]
            for metricIndex in xrange(0,metricsCount):
                returnMetricsMatrixRet[metricIndex].append(Metric.funcGetAlphaMetric(ldAbundancies = sampleAbundance, strMetric = lsDiversityMetricAlpha[metricIndex]))
        return returnMetricsMatrixRet

    #Testing 6 cases
    @staticmethod
    def funcGetBetaMetric(npadAbundancies=None, sMetric=None, istrmTree=None, istrmEnvr=None, lsSampleOrder=None):
        """
        Takes a matrix of values and returns a beta metric matrix. The metric returned is indicated by name (sMetric).
		
        :param	npadAbundancies:	Numpy array of sample abundances to measure against.
        :type:	Numpy Array	Numpy array where row=samples and columns = features.
        :param	sMetric:	String name of beta metric. Possibilities are listed in microPITA.
        :type:	String	String name of beta metric. Possibilities are listed in microPITA.
        :return	Double:	Measurement indicated by metric for given abundance list
        """
	
        if sMetric == Metric.c_strBrayCurtisDissimilarity:
            return Metric.funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies=npadAbundancies)
        elif sMetric == Metric.c_strInvBrayCurtisDissimilarity:
            return Metric.funcGetInverseBrayCurtisDissimilarity(ldSampleTaxaAbundancies=npadAbundancies)
        elif sMetric in Metric.setBetaDiversities:
            return Metric.funcGetDissimilarityByName(ldSampleTaxaAbundancies=npadAbundancies, strMetric=sMetric)
        elif sMetric == Metric.c_strUnifracUnweighted:
            xReturn = Metric.funcGetUnifracDistance(istrmTree=istrmTree,istrmEnvr=istrmEnvr,lsSampleOrder=lsSampleOrder,fWeighted=False)
            return xReturn[0] if not False else xReturn
        elif sMetric == Metric.c_strUnifracWeighted:
            xReturn = Metric.funcGetUnifracDistance(istrmTree=istrmTree,istrmEnvr=istrmEnvr,lsSampleOrder=lsSampleOrder,fWeighted=True)
            return xReturn[0] if not False else xReturn
        else:
            return False

    #Test Cases 11
    @staticmethod
    def funcReadMatrixFile(istmMatrixFile, lsSampleOrder=None):
	"""
	Reads in a file with a precalculated beta-diversty matrix.

	:param istmMatrixFile:	File with beta-diversity matrix
	:type:	FileStream of String file path
	"""

        #Read in data
        f = csv.reader(open(istmMatrixFile,"r") if isinstance(istmMatrixFile, str) else istmMatrixFile, delimiter=ConstantsBreadCrumbs.c_matrixFileDelim )

        #Get header
        try:
            lsHeader = f.next()
        except StopIteration:
            return ([],[])
        lsHeaderReducedToSamples = [sHeader for sHeader in lsHeader if sHeader in lsSampleOrder] if lsSampleOrder else lsHeader[1:]

        #If no sample ordering is given, set the ordering to what is in the file
        if not lsSampleOrder:
	    lsSampleOrder = lsHeaderReducedToSamples

        #Preallocate matrix
        mtrxData = np.zeros(shape=(len(lsSampleOrder),len(lsSampleOrder)))

        #Make sure all samples requested are in the file
        if(not len(lsSampleOrder) == len(lsHeaderReducedToSamples)): return False

	for lsLine in f:
            if lsLine[0] in lsSampleOrder:
                iRowIndex = lsSampleOrder.index(lsLine[0])

                for i in xrange(1,len(lsSampleOrder)):
                    iColumnIndexComing = lsHeader.index(lsSampleOrder[i])
                    iColumnIndexGoing = lsSampleOrder.index(lsSampleOrder[i])
                    mtrxData[iRowIndex,iColumnIndexGoing] = lsLine[iColumnIndexComing]
                    mtrxData[iColumnIndexGoing,iRowIndex] = lsLine[iColumnIndexComing]
        tpleMData = mtrxData.shape
        mtrxData = mtrxData if any(sum(ld)>0 for ld in mtrxData) or ((tpleMData[0]==1) and (tpleMData[1]==1)) else []
        return (mtrxData,lsSampleOrder)

    #Test cases 2
    @staticmethod
    def funcWriteMatrixFile(mtrxMatrix, ostmMatrixFile, lsSampleNames=None):
        """
        Writes a square matrix to file.
        
        :param mtrxMatrix:	Matrix to write to file
        :type:	Numpy array
        :lsSampleNames:	The names of the samples in the order of the matrix
        :type:	List of strings
        :ostmBetaMatrixFile:	File to write to
        :type:	String or file stream
        """

        if not sum(mtrxMatrix.shape)>0 or not ostmMatrixFile:
            return False

        #Check to make sure the sample names are the correct length
        tpleiShape = mtrxMatrix.shape
        if not lsSampleNames:
            lsSampleNames = range(tpleiShape[0])
        if not(len(lsSampleNames) == tpleiShape[0]):
            print "".join(["Metric.funcWriteMatrixFile. Error= Length of sample names ("+str(len(lsSampleNames))+") and matrix ("+str(mtrxMatrix.shape)+") not equal."])
            return False

        #Write to file
        ostmOut = csv.writer(open(ostmMatrixFile,"w") if isinstance(ostmMatrixFile,str) else ostmMatrixFile, delimiter=ConstantsBreadCrumbs.c_matrixFileDelim )

        #Add the additional space at the beginning of the sample names to represent the id row/column
        lsSampleNames = [""]+lsSampleNames

        #Write header and each row to file
        ostmOut.writerow(lsSampleNames)
        [ostmOut.writerow([lsSampleNames[iIndex+1]]+mtrxMatrix[iIndex,].tolist()) for iIndex in xrange(tpleiShape[0])]
        return True
