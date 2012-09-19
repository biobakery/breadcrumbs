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
import numpy as np
from types import *
from ValidateData import ValidateData

#External libraries
from cogent.maths.stats.alpha_diversity import chao1_uncorrected, chao1_bias_corrected
from scipy.spatial.distance import pdist

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

    #Additive inverses of beta metrics
    c_strInvBrayCurtisDissimilarity = "InB_Curtis"

    #Richness
    c_strShannonRichness = "ShannonR"
    c_strObservedCount = "Observed_Count"

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
            return chao1_bias_corrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)
        else:
            return chao1_uncorrected(observed = totalObservedSpecies, singles = singlesObserved, doubles = doublesObserved)

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
            return pdist(ldSampleTaxaAbundancies, funcDistanceFunction)
        except ValueError as error:
            print "".join(["Metric.funcGetDissimilarity. Error=",str(error)])
            return False

    #Test 3
    @staticmethod
    def funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies=None):
        """
        Calculates the BrayCurtis Beta dissimilarity index.
        d(u,v)=sum(abs(row1-row2))/sum(row1+row2).
        This is scale invariant.
        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:
        :type:	List of doubles
        :return	Double:	Dissimilarity metric
        """

        #Calculate metric
        try:
            return pdist(X=ldSampleTaxaAbundancies, metric='braycurtis')
        except ValueError as error:
            print "".join(["Metric.getBrayCurtisDissimilarity. Error=",str(error)])
            return False

    #Test 3
    @staticmethod
    def funcGetInverseBrayCurtisDissimilarity(ldSampleTaxaAbundancies=None):
        """
        Calculates 1 - the BrayCurtis Beta dissimilarity index.
        d(u,v)=1-(sum(abs(row1-row2))/sum(row1+row2)).
        This is scale invariant and ranges between 0 and 1.
        If you have 5 rows (labeled r1,r2,r3,r4,r5) the vector are the distances in this order.
        condensed form = [d(r1,r2), d(r1,r3), d(r1,r4), d(r1,r5), d(r2,r3), d(r2,r4), d(r2,r5), d(r3,r4), d(r3,r5), d(r4,r5)].
        Note***: Assumes that the abundance measurements are already normalized by the total population N.

        :param	ldSampleTaxaAbundancies:	An np.array of samples (rows) x measurements (columns) in which distance is measured between rows
        :type:	List	List of doubles
        :return	Double	1 - Bray-Curtis dissimilarity.	
        """

        bcValue = Metric.funcGetBrayCurtisDissimilarity(ldSampleTaxaAbundancies = ldSampleTaxaAbundancies)
        if not type(bcValue) is BooleanType:
            return 1.0-bcValue
        return False

    #Test 4
    @staticmethod
    def funcGetAlphaMetric(ldAbundancies=None, strMetric=None):
        """
        Get alpha abundance of the metric for the vector.

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
            return Metric.funcGetObservedCount(npaSampleAbundances=ldAbundancies)
        #Needs NOT Normalized Abundance
        elif(strMetric == Metric.c_strChao1Diversity):
            return Metric.funcGetChao1DiversityIndex(ldSampleTaxaAbundancies=ldAbundancies)
        else:
            return False

    #Test 3
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
