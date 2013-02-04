"""
Author: Timothy Tickle
Description: Class to abstract an abundance table and methods to run on such a table.
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

import csv
import sys
import blist
from CClade import CClade
from ConstantsBreadCrumbs import ConstantsBreadCrumbs
import copy
import numpy as np
import os
import re
import scipy.stats
import string
from ValidateData import ValidateData

c_dTarget	= 1.0
c_fRound	= False
c_iSumAllCladeLevels = -1
c_fOutputLeavesOnly = False

class AbundanceTable:
	"""
	Represents an abundance table and contains common function to perform on the object.

	This class is made from an abundance data file. What is expected is a text file delimited by
	a character (which is given to the object). The first column is expected to be the id column
	for each of the rows. Metadata is expected before measurement data. Columns are samples and
	rows are features (bugs). 
	"""

	def __init__(self, npaAbundance, dictMetadata, strName, strLastMetadata, lOccurenceFilter = None, cFileDelimiter = ConstantsBreadCrumbs.c_cTab, cFeatureNameDelimiter="|"):
		"""
		Constructor for an abundance table.

		:param	npaAbundance:	Structured Array of abundance data (Row=Features, Columns=Samples)
		:type:	Numpy Structured Array abundance data (Row=Features, Columns=Samples)
		:param	dictMetadata:	Structured Array of abundance data (Row=Features, Columns=Samples)
		:type:	Dictionary	Dictionary of metadata {"String ID":["strValue","strValue","strValue","strValue","strValue"]}
	 	:param	strName:	The name of the metadata that serves as the ID for the columns (For example a sample ID)
		:type:	string
		:param	strLastMetadata: The string last metadata name
      		:type:	string
		:param	lOccurenceFilter: List of integers used in an occurence filter. [Min abundance, Min sample]
		:type:	List of integers
		:param	cFileDelimiter:	Character used as the delimiter of the file that is read in to create the abundance table.
								Will also be used to write the abudance table file to a file to keep file consistency.
		:type:	Character delimiter for reading the data in (default = TAB)
		:param	cFeatureNameDelimiter:	Character used as the delimiter of the feature names (column 1). This is useful if the name are complex, for instance consensus lineages in metagenomics.
		:type:	Character delimiter for feature names (default = |)
		"""

		#The abundance data
		self._npaFeatureAbundance = npaAbundance

		#The metdata
		self._dictTableMetadata = dictMetadata

		#The name of the object relating to the file it was read from or would have been read from if it exists
		#Keeps tract of changes to the file through the name
		#Will be used to wrtie out the object to a file as needed
		self._strOriginalName = strName

		#The original number of features in the table
		self._iOriginalFeatureCount = -1

		#The original number of samples in the table
		self._iOriginalSampleCount = -1

		#Indicates if the table has been filtered and how
		self._strCurrentFilterState = ""

		#The feature name delimiter
		self._cFeatureDelimiter = cFeatureNameDelimiter

		#The delimiter from the source file
		self._cDelimiter = cFileDelimiter

		#The lastmetadata name (which should be preserved when writing the file)
		self._strLastMetadataName = strLastMetadata

		self._fIsNormalized = self._fIsSummed = None
		#If contents is not a false then set contents to appropriate objects
		if ( self._npaFeatureAbundance != None ) and self._dictTableMetadata:
			self._iOriginalFeatureCount = self._npaFeatureAbundance.shape[0]
			self._iOriginalSampleCount = len(self.funcGetSampleNames())
		
			self._fIsNormalized = ( max( [max( list(a)[1:] or [0] ) for a in self._npaFeatureAbundance] or [0] ) <= 1 )
			lsLeaves = AbundanceTable.funcGetTerminalNodesFromList( [a[0] for a in self._npaFeatureAbundance], self._cFeatureDelimiter )
			self._fIsSummed = ( len( lsLeaves ) != len( self._npaFeatureAbundance ) )

			#Occurence filtering
			#Removes features that do not have a given level iLowestAbundance in a given amount of samples iLowestSampleOccurence
			if ( not self._fIsNormalized ) and lOccurenceFilter:
				iLowestAbundance, iLowestSampleOccurrence = lOccurenceFilter
				self.funcFilterAbundanceBySequenceOccurence( iLowestAbundance, iLowestSampleOccurrence )
#	  else:
#		sys.stderr.write( "Abundance or metadata was None, should be atleast an empty object\n" )

	@staticmethod
	def funcMakeFromFile(xInputFile, cDelimiter = ConstantsBreadCrumbs.c_cTab, sMetadataID = None, sLastMetadata = None,
	   lOccurenceFilter = None, cFeatureNameDelimiter="|", xOutputFile = None):
		"""
		Creates an abundance table from a table file.

		:param	xInputFile:	Path to input file.
		:type:	String		String path.
		:param	cDelimiter:	Delimiter for parsing the input file.
		:type:	Character	Character.
		:param	sMetadataID:	String ID that is a metadata row ID (found on the first column) and used as an ID for samples
		:type:	String		String ID
		:param	sLastMetadata:	The ID of the metadata that is the last metadata before measurement or feature rows.
		:type:	String		String ID
		:param	lOccurenceFilter: List of integers used in an occurence filter. [Min abundance, Min sample]
		:type:	List of integers
		:param	cFeatureNameDelimiter:	Used to parse feature (bug) names if they are complex.
						For example if they are consensus lineages and contain parent clade information.
		:type:	Character	Delimiting letter
		:param	xOutputFile:	File to output the abundance table which was read in.
		:type:	FileStream or String file path
		:return	AbundanceTable:	Will return an AbundanceTable object on no error. Returns False on error.
		"""
		
		#Get output file and remove if existing
		outputFile = open( strOutputFileName, "w" ) if isinstance(xOutputFile, str) else xOutputFile
		#Read in from text file to create the abundance and metadata structures
		lContents = AbundanceTable._funcTextToStructuredArray(xInputFile=xInputFile, cDelimiter=cDelimiter,
				sMetadataID = sMetadataID, sLastMetadata = sLastMetadata, ostmOutputFile = outputFile)

		#If contents is not a false then set contents to appropriate objects
		return AbundanceTable(npaAbundance=lContents[0], dictMetadata=lContents[1], strName=str(xInputFile), strLastMetadata=sLastMetadata,
		  lOccurenceFilter = lOccurenceFilter, cFileDelimiter=cDelimiter, cFeatureNameDelimiter=cFeatureNameDelimiter) if lContents else False

	#Testing Status: Light happy path testing
	@staticmethod
	def funcCheckRawDataFile(strReadDataFileName, iFirstDataIndex = -1, sLastMetadataName = None, lOccurenceFilter = None, strOutputFileName = "", cDelimiter = ConstantsBreadCrumbs.c_cTab):
		"""
		Check the input abundance table.
		Currently reduces the features that have no occurence.
		Also inserts a NA for blank metadata and a 0 for blank abundance data.
		Gives the option to filter features through an occurence filter (a feature must have a level of abundance in a minimal number of samples to be included).
		Either iFirstDataIndex or sLastMetadataName must be given

		:param	strReadDataFileName:	File path of file to read and check.
		:type:	String	File path.
		:param	iFirstDataIndex:	First (row) index of data not metadata in the abundance file.
		:type:	Integer	Index starting at 0.
		:param	sLastMetadataName:	The ID of the last metadata in the file. Rows of measurements should follow this metadata.
		:type:	String
		:param	lOccurenceFilter:	The lowest number of occurences in the lowest number of samples needed for a feature to be kept
		:type:	List[2]	List length 2 [lowest abundance (not normalized), lowest number of samples to occur in] (eg. [2.0,2.0])
		:param	strOutputFileName:	File path of out put file.
		:type:	String	File path.
		:param	cDelimiter:	Character delimiter for reading and writing files.
		:type:	Character	Delimiter.
		:return	Output Path:	Output path for written checked file.
		"""

		#Validate parameters
		if (iFirstDataIndex == -1) and (sLastMetadataName == None):
			print "AbundanceTable:checkRawDataFile::Error, either iFirstDataIndex or sLastMetadataNamemust be given."
			return False

		#Get output file and remove if existing
		outputFile = strOutputFileName
		if not strOutputFileName:
			outputFile = os.path.splitext(strReadDataFileName)[0]+ConstantsBreadCrumbs.OUTPUT_SUFFIX

		#Read input file lines
		#Drop blank lines
		readData = ""
		with open(strReadDataFileName,'r') as f:
			readData = f.read()
		readData = filter(None,readData.split(ConstantsBreadCrumbs.c_strEndline))

		#Read the length of each line and make sure there is no jagged data
		#Also hold row count for the metadata
		iLongestLength = len(readData[0].split(cDelimiter))
		iMetadataRow = -1
		if not sLastMetadataName:
			sLastMetadataName = "None"
		for iIndex, strLine in enumerate(readData):
			sLineElements = strLine.split(cDelimiter)
		if sLineElements[0] == sLastMetadataName:
			iMetadataRow = iIndex
		iLongestLength = max(iLongestLength, len(sLineElements))

		#If not already set, set iFirstDataIndex
		if iFirstDataIndex < 0:
			iFirstDataIndex = iMetadataRow + 1

		#Used to substitute . to -
		reSubPeriod = re.compile('\.')

		#File writer
		with open(outputFile,'w') as f:

			#Write metadata
			#Empty data is changed to a default
			#Jagged ends are filled with a default
			for strDataLine in readData[:iFirstDataIndex]:
				lsLineElements = strDataLine.split(cDelimiter)
				for iindex, sElement in enumerate(lsLineElements):
					if not sElement.strip():
						lsLineElements[iindex] = ConstantsBreadCrumbs.c_strEmptyDataMetadata
				if len(lsLineElements) < iLongestLength:
					lsLineElements = lsLineElements + ([ConstantsBreadCrumbs.c_strEmptyDataMetadata]*(iLongestLength-len(lsLineElements)))
				f.write(cDelimiter.join(lsLineElements)+ConstantsBreadCrumbs.c_strEndline)

			#For each data line in the table
			for line in readData[iFirstDataIndex:]:
				writeToFile = False
				cleanLine = list()
				#Break line into delimited elements
				lineElements = line.split(cDelimiter)

				#Clean feature name
				sCleanFeatureName = reSubPeriod.sub("-",lineElements[0])

				#For each element but the first (taxa name)
				#Element check to see if not == zero
				#If so add to output
				for element in lineElements[1:]:
					if(element.strip() in string.whitespace):
						cleanLine.append(ConstantsBreadCrumbs.c_strEmptyAbundanceData)
					#Set abundance of 0 but do not indicate the line should be saved
					elif(element == "0"):
						cleanLine.append(element)
					#If an abundance is found set the line to be saved.
					else:
						cleanLine.append(element)
						writeToFile = True

				#Occurence filtering
				#Removes features that do not have a given level iLowestAbundance in a given amount of samples iLowestSampleOccurence
				if lOccurenceFilter:
					iLowestAbundance, iLowestSampleOccurence = lOccurenceFilter
					if iLowestSampleOccurence > sum([1 if float(sEntry) >= iLowestAbundance else 0 for sEntry in cleanLine]):
						writeToFile = False

				#Write to file
				if writeToFile:    
					f.write(sCleanFeatureName+cDelimiter+cDelimiter.join(cleanLine)+ConstantsBreadCrumbs.c_strEndline)
		return outputFile

	def __repr__(self):
		"""
		Represent or print object.
		"""
		return "AbundanceTable"

	def __str__(self):
	  """
	  Create a string representation of the Abundance Table.
	  """

	  return "".join(["Sample count:", str(len(self._npaFeatureAbundance.dtype.names[1:])),
	  os.linesep+"Feature count:", str(len(self._npaFeatureAbundance[self._npaFeatureAbundance.dtype.names[0]])),
	  os.linesep+"Id Metadata:", self._npaFeatureAbundance.dtype.names[0],
	  os.linesep+"Metadata ids:", str(self._dictTableMetadata.keys()),
	  os.linesep+"Metadata count:", str(len(self._dictTableMetadata.keys())),
	  os.linesep+"Originating source:",self._strOriginalName,
	  os.linesep+"Original feature count:", str(self._iOriginalFeatureCount),
	  os.linesep+"Original sample count:", str(self._iOriginalSampleCount),
	  os.linesep+"Is normalized:", str(self._fIsNormalized),
	  os.linesep+"Is summed:", str(self._fIsSummed),
	  os.linesep+"Current filtering state:", str(self._strCurrentFilterState),
	  os.linesep+"Feature delimiter:", self._cFeatureDelimiter,
	  os.linesep+"File delimiter:",self._cDelimiter])

	#Testing Status: Light happy path testing
	@staticmethod
	def _funcTextToStructuredArray(xInputFile = None, cDelimiter = ConstantsBreadCrumbs.c_cTab, sMetadataID = None, sLastMetadata = None, ostmOutputFile = None):
		"""
		Private method
		Used to read in a file that is samples (column) and taxa (rows) into a structured array.

		:param	xInputFile:	File stream or path to input file.
		:type:	String		File stream or string path.
		:param	cDelimiter:	Delimiter for parsing the input file.
		:type:	Character	Character.
		:param	sMetadataID:	String ID that is a metadata row ID (found on the first column) and used as an ID for samples
		:type:	String		String ID
		:param	sLastMetadata:	The ID of the metadata that is the last metadata before measurement or feature rows.
		:type:	String		String ID
		:param	ostmOutputFile:	Output File to write to if needed. None does not wrtie the file.
		:type:	FileStream or String
		:return	[taxData,metadata]:	Numpy Structured Array of abundance data and dictionary of metadata.
										Metadata is a dictionary as such {"ID", [value,value,values...]}
										Values are in the order thety are read in (and the order of the sample names).
										ID is the first column in each metadata row.
										[Numpy structured Array, Dictionary]
		"""

		istmInput = open( xInputFile ) if isinstance(xInputFile, str) else xInputFile
		iFirstDataRow = -1
		namesRow = None
		metadata = dict()
		dataMatrix = []
		iIndex = -1
		csvw = None
		if ostmOutputFile:
			csvw = csv.writer( open(ostmOutputFile,'w') if isinstance(ostmOutputFile, str) else ostmOutputFile, csv.excel_tab, delimiter = cDelimiter )
		for lsLineElements in csv.reader( istmInput, csv.excel_tab, delimiter = cDelimiter ):
			iIndex += 1
			taxId, sampleReads = lsLineElements[0], lsLineElements[1:]
			#Data
			if iFirstDataRow > 0:
				try:
					dataMatrix.append(tuple([taxId] + [( float(s) if s.strip( ) else 0 ) for s in sampleReads]))
				except ValueError:
					sys.stderr.write( "AbundanceTable:textToStructuredArray::Error, non-numerical value on data row. File:" + str(xInputFile) +
						" Row:" + str(lsLineElements) + "\n" )
					return False
			#Metadata
			else:
				for i, s in enumerate( sampleReads ):
					if not s.strip( ):
						sampleReads[i] = ConstantsBreadCrumbs.c_strEmptyDataMetadata
				if ( ( not sMetadataID ) and ( iIndex == 0 ) ) or ( taxId == sMetadataID ):
					namesRow = lsLineElements
				metadata[taxId]=sampleReads
				if ( not sLastMetadata ) or ( taxId == sLastMetadata ):
					iFirstDataRow = iIndex + 1
			if csvw:
				csvw.writerow( [taxId] + sampleReads )
		if sLastMetadata and ( not dataMatrix ):
			sys.stderr.write( "AbundanceTable:textToStructuredArray::Error, did not find the row for the last metadata ID. File:" + str(xInputFile) +
				" Identifier:" + sLastMetadata + "\n" )
			return False
		#Make sure the names are found
		if namesRow == None:
			sys.stderr.write( "AbundanceTable:textToStructuredArray::Error, did not find the row for the unique sample/column. File:" + str(xInputFile) +
				" Identifier:" + sMetadataID + "\n" )
			return False

		#Now we know the longest taxId we can define the first column holding the tax id
		#Gross requirement of Numpy structured arrays, a = ASCII followed by max # of characters (as a string)
		longestTaxId = max( len(a[0]) for a in dataMatrix )
		dataTypeVector = [(namesRow[0],'a' + str(longestTaxId*2))] + [(s, "f4") for s in namesRow[1:]]
		#Create structured array
		taxData = np.array(dataMatrix,dtype=np.dtype(dataTypeVector))

		return [taxData,metadata]

	#Happy path tested
	def funcGetSampleNames(self):
		"""
		Returns the sample names (IDs) contained in the abundance table.

		:return	Sample Name:	A List of sample names indicated by the metadata associated with the sMetadataId given in table creation.
								A list of string names or empty list on error as well as no underlying table.
		"""

		return self._npaFeatureAbundance.dtype.names[1:] if ( self._npaFeatureAbundance != None ) else []

	#Happy Path Tested
	def funcGetIDMetadataName(self):
		"""
		Returns the metadata id.

		:return	ID:	The metadata id (the sample Id).
					  Returns none on error.
		"""

		return self._npaFeatureAbundance.dtype.names[0] if ( self._npaFeatureAbundance != None ) else None

	#Happy path tested
	def funcGetAbundanceCopy(self):
		"""
		Returns a deep copy of the abundance table.

		:return	Numpy Structured Array:	The measurement data in the Abundance table. Can use sample names to access each column of measurements.
									   Returns none on error.
		"""

		return self._npaFeatureAbundance.copy() if ( self._npaFeatureAbundance != None ) else None

	#Happy path tested
	def funcGetAverageAbundancePerSample(self, lsTargetedFeatures):
		"""
		Averages feature abundance within a sample.
	
		:param	lsTargetedFeatures:	String names of features to average
		:type:	List of string names of features which are measured
		:return	List: List of lists or boolean (False on error). One internal list per sample indicating the sample and the feature's average abudance
			[[sample,average abundance of selected taxa]] or False on error
		"""

		#Sample rank averages [[sample,average abundance of selected taxa]]
		sampleAbundanceAverages = []
		
		sampleNames = self.funcGetSampleNames()
		allTaxaNames = self.funcGetFeatureNames()
		#Get an abundance table compressed to features of interest
		abndReducedTable = self.funcGetFeatureAbundanceTable(lsTargetedFeatures)
		if abndReducedTable == None:
			return False  

		#If the taxa to be selected are not in the list, Return nothing and log
		lsMissing = []
		for sFeature in lsTargetedFeatures:
			if not sFeature in allTaxaNames:
				lsMissing.append(sFeature)
			else:
				#Check to make sure the taxa of interest is not average abundance of 0
				if not abndReducedTable.funcGetFeatureSumAcrossSamples(sFeature):
					lsMissing.append(sFeature)
		if len(lsMissing) > 0:
			sys.stderr.write( "Could not find features for averaging: " + str(lsMissing) )
			return False

		#For each sample name get average abundance
		for sName in sampleNames:
			npaFeaturesSample = abndReducedTable.funcGetSample(sName)
			sampleAbundanceAverages.append([sName,sum(npaFeaturesSample)/float(len(npaFeaturesSample))])

		#Sort based on average
		return sorted(sampleAbundanceAverages, key = lambda sampleData: sampleData[1], reverse = True)

	#Happy path tested 1
	def funcGetAverageSample(self):
		"""
		Returns the average sample of the abundance table.
		This average sample is made of the average of each feature.
		:return list: A list of averages in the order of the feature names.
		"""

		ldAverageSample = []
		#If there are no samples then return empty list.
		if len(self.funcGetSampleNames()) < 1:
			return ldAverageSample

		#If there are samples return the average of each feature in the order of the feature names.
		for sFeature in self._npaFeatureAbundance:
			npFeaturesAbundance = list(sFeature)[1:]
			ldAverageSample.append(sum(npFeaturesAbundance)/float(len(npFeaturesAbundance)))

		return ldAverageSample

	#Happy Path Tested
	def funcGetFeatureAbundanceTable(self, lsFeatures):
		"""
		Returns a copy of the current abundance table with the abundance of just the given features.

		:param	lsFeatures:	String Feature IDs that are kept in the compressed abundance table.
		:type:	List of strings	Feature IDs (found as the first entry of a filter in the input file.
		:return	AbundanceTable:	A compressed version of the abundance table.
				  On an error None is returned.
		"""
		
		if ( self._npaFeatureAbundance == None ) or ( lsFeatures == None ):
			return None

		#Get a list of boolean indicators that the row is from the features list
		lfFeatureData = [sRowID in lsFeatures for sRowID in self.funcGetFeatureNames()]
		#compressed version as an Abundance table
		lsNamePieces = os.path.splitext(self._strOriginalName)
		return AbundanceTable(npaAbundance=np.compress(lfFeatureData, self._npaFeatureAbundance, axis = 0),
					dictMetadata = self.funcGetMetadataCopy(),
					strName = lsNamePieces[0] + "-" + str(len(lsFeatures)) +"-Features"+lsNamePieces[1],
					strLastMetadata=self.funcGetLastMetadataName(),
					cFileDelimiter = self.funcGetFileDelimiter(), cFeatureNameDelimiter= self.funcGetFeatureDelimiter())

	#Happy path tested
	def funcGetFeatureDelimiter(self):
		"""
		The delimiter of the feature names (For example to use on concensus lineages).

		:return	Character:	Delimiter for the feature name pieces if it is complex.
		"""

		return self._cFeatureDelimiter

	#Happy path tested
	def funcGetFeatureCount(self):
		"""
		Returns the current feature count.

		:return	Count:	Returns the int count of features in the abundance table.
						Returns None on error.
		"""

		return self._npaFeatureAbundance.shape[0] if not self._npaFeatureAbundance is None else 0

	#Happy path tested
	def funcGetFeatureSumAcrossSamples(self,sFeatureName):
		"""
		Returns float sum of feature values across the samples.

		:param	sFeatureName: The feature ID to get the sum from.
		:type:	String.
		:return	Double:	Sum of one feature across samples.
		"""

		for sFeature in self._npaFeatureAbundance:
			if sFeature[0] == sFeatureName:
				return sum(list(sFeature)[1:])
		return None

	#Happy path tested
	def funcGetFeatureNames(self):
		"""
		Return the feature names as a list.

		:return	Feature Names:	List of feature names (or IDs) as strings.
								As an error returns empty list.
		"""

		if (not self._npaFeatureAbundance == None):
			return self._npaFeatureAbundance[self.funcGetIDMetadataName()]
		return []

	#Happy path tested
	def funcGetFileDelimiter(self):
		"""
		The delimiter of the file the data was read from and which is also the delimiter which would be used to write the data to a file.

		:return	Character:	Delimiter for the parsing and writing the file.
		"""

		return self._cDelimiter

	def funcGetLastMetadataName(self):
		"""
		Get the last metadata name that seperates abundance and metadata measurements.

		:return string:	Metadata name
		"""
		return self._strLastMetadataName

	#Happy path tested
	def funcGetSample(self,sSampleName):
		"""
		Return a copy of the feature measurements of a sample.

		:param	sSampleName:	Name of sample to return.	
		:type:	String	
		:return	Sample: Measurements	Feature measurements of a sample.
				Empty numpy array returned on error.
		"""

		if (not self._npaFeatureAbundance == None):
			return self._npaFeatureAbundance[sSampleName].copy()
		return np.array([])

	#Happy path tested
	def funcGetMetadata(self, strMetadataName):
		"""
		Returns a list of metadata that is associated with the given metadata name (id).

		:param	strMetadataName:	String metadata ID to be returned
		:type:	String	ID
		:return	Metadata:	List of metadata
		"""
		
		return copy.deepcopy( self._dictTableMetadata.get(strMetadataName) ) \
			if self._dictTableMetadata else None

	#Happy path tested
	def funcGetMetadataCopy(self):
		"""
		Returns a deep copy of the metadata.

		:return	Metadata copy:	{"ID":[value,value...]}
		"""

		return copy.deepcopy(self._dictTableMetadata)
		
	#Happy path tested
	def funcGetName(self):
		"""
		Returns the name of the object which is the file name that generated it.
		If the object was generated from an Abundance Table (for instance through stratification)
		the name is still in the form of a file that could be written to which is informative
		of the changes that have occurred on the data set.
		:return string: Name
		"""
		return self._strOriginalName

	#Happy path tested. could do more
	def funcGetTerminalNodes(self):
		"""
		Returns the terminal nodes given the current feature names in the abundance table. The 
		features must contain a consensus lineage or all will be returned.
		:return List:	List of strings of the terminal nodes given the abundance table.
		"""
		return AbundanceTable.funcGetTerminalNodesFromList(lsNames=self.funcGetFeatureNames(),cNameDelimiter=self.funcGetFeatureDelimiter())

	#Tested 2 test cases
	@staticmethod
	def funcGetTerminalNodesFromList(lsNames,cNameDelimiter):
		"""
		Returns the terminal nodes given the current feature names in the abundance table. The 
		features must contain a consensus lineage or all will be returned.

		:param	lsNames:	The list of string names to parse and filter.
		:type:	List of strings
		:param	cNameDelimiter:	The delimiter for the name of the features.
		:type:	Character	Delimiter
		:return list:	A list of terminal elements in the list (given only the list).
		"""

		#Build hash
		dictCounts = dict()
		for strTaxaName in lsNames:
			#Split into the elements of the clades
			lsClades = filter(None,strTaxaName.split(cNameDelimiter))
			#Count clade levels
			iCladeLength = len(lsClades)

			#Evaluate first element
			sClade = lsClades[0]
			dictCounts[sClade] = sClade not in dictCounts

			#Evaluate the rest of the elements
			if iCladeLength < 2:
				continue
			for iIndex in xrange(1,iCladeLength):
				prevClade = sClade
				sClade = cNameDelimiter.join([sClade,lsClades[iIndex]])
				if sClade in dictCounts:
					dictCounts[sClade] = dictCounts[prevClade] = False
				else:
					dictCounts[sClade] = True
					dictCounts[prevClade] = False

		#Return only the elements that were of count 1
		return filter( lambda s: dictCounts[s] == True, dictCounts )

	#Happy path tested
	def funcIsNormalized(self):
		"""
		Returns if the data has been normalized.

		:return	Boolean:	Indicates if the data is normalized.
						   True indicates it the data is normalized.
		"""

		return self._fIsNormalized

	#Happy path tested
	def funcIsPrimaryIdMetadata(self,sMetadataName):
		"""
		Checks the metadata data associatd with the sMetadatName and returns if the metadata is unique.
		This is important to some of the functions in the Abundance Table specifically when translating from one metadata to another.
		
		:param	sMetadataName:	ID of metadata to check for uniqueness.
		:type:	String	Metadata ID.
		:return	Boolean:	Returns indicator of uniqueness.
							True indicates unique.
		"""

		lMetadata = self.funcGetMetadata(sMetadataName)
		if not lMetadata:
			return False
		return (len(lMetadata) == len(set(lMetadata)))

	#Happy path tested
	def funcIsSummed(self):
		"""
		Return is the data is summed.

		:return	Boolean:	Indicator of being summed. True indicates summed.
		"""

		return self._fIsSummed

	#Happy path tested
	def funcFilterAbundanceByPercentile(self, dPercentileCutOff = 95.0, dPercentageAbovePercentile=1.0):
		"""
		Filter on features.
		A feature is removed if it's abundance is not found in the top X percentile a certain percentage of the samples.

		:param	dPercentileCutOff:	The percentile used for filtering.
		:type:	double	A double between 0.0 and 100.0
		:param	dPercentageAbovePercentile:	The percentage above the given percentile (dPercentileCutOff) that must exist to keep the feature.
		:type:	double	Between 0.0 and 100.0
		:return	Boolean:	Indicator of filtering occuring without error. True indicates filtering occuring.
		"""

		#No need to do anything
		if(dPercentileCutOff==0.0) or (dPercentageAbovePercentile==0.0):
			return True

		#Sample names
		lsSampleNames = self.funcGetSampleNames()

		#Scale percentage out of 100
		dPercentageAbovePercentile = dPercentageAbovePercentile/100.0

		#Sample count
		iSampleCount = len(lsSampleNames)

		#Get a threshold score of the value at the specified percentile for each sample
		#In the order of the sample names
		ldScoreAtPercentile = [scipy.stats.scoreatpercentile(self._npaFeatureAbundance[lsSampleNames[iIndex]],dPercentileCutOff) for iIndex in xrange(iSampleCount)]

		#Record how many entries for each feature have a value equal to or greater than the dPercentileCutOff
		#If the percentile of entries passing the criteria are above the dPercentageAbovePercentile put index in list to keep
		liKeepIndices = []
		iSampleCount = float(iSampleCount)
		for iRowIndex, npaRow in enumerate(self._npaFeatureAbundance):
			iCountPass = sum([1 if dValue >= ldScoreAtPercentile[iValueIndex] else 0 for iValueIndex, dValue in enumerate(list(npaRow)[1:])])
			if (iCountPass / iSampleCount) >= dPercentageAbovePercentile:
				liKeepIndices.append(iRowIndex)

		#Compress array
		self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepIndices,:]

		#Update filter state
		self._strCurrentFilterState += ":dPercentileCutOff=" + str(dPercentileCutOff) + ",dPercentageAbovePercentile=" + str(dPercentageAbovePercentile)

		return True

	#Happy path tested
	def funcFilterAbundanceBySequenceOccurence(self, iMinSequence = 2, iMinSamples = 2):
		"""
		Filter abundance by requiring features to have a minimum sequence occurence in a minimum number of samples.
		Will evaluate greater than or equal to the iMinSequence and iMinSamples.

		:param	iMinSequence:	Minimum sequence to occur.
		:type:	Integer	Number Greater than 1.
		:param	iMinSamples:	Minimum samples to occur in.
		:type:	Integer	Number greater than 1.
		:return	Boolean:	Indicator of the filter running without error. False indicates error.
		"""

		#No need to do anything
		if(iMinSequence==0) or (iMinSamples==0):
			return True

		#This normalization requires the data to be reads
		if self._fIsNormalized:
			#sys.stderr.write( "Could not filter by sequence occurence because the data is already normalized.\n" )
			return False

		#Holds which indexes are kept
		liKeepFeatures = []
		for iRowIndex, dataRow in enumerate( self._npaFeatureAbundance ):
			#See which rows meet the criteria and keep the index if needed.
			if len( filter( lambda d: d >= iMinSequence, list(dataRow)[1:] ) ) >= iMinSamples:
				liKeepFeatures.append(iRowIndex)

		#Compress array
		self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]
		#Update filter state
		self._strCurrentFilterState += ":iMinSequence=" + str(iMinSequence) + ",iMinSamples=" + str(iMinSamples)

		return True
   
	#1 Happy path test
	def funcFilterFeatureBySD(self, dMinSDCuttOff = 0.0):
		"""
		A feature is removed if it's abundance is not found to have standard deviation more than the given dMinSDCutoff.

		:param	dMinSDCuttOff:	Standard deviation threshold.
		:type:	Double	A double greater than 0.0.
		:return	Boolean:	Indicator of success. False indicates error.
		"""

		#No need to do anything
		if(dMinSDCuttOff==0.0):
			return True

		#Holds which indexes are kept
		liKeepFeatures = []

		#Evaluate each sample
		for iRowIndex, dataRow in enumerate(self._npaFeatureAbundance):
			if(np.std(list(dataRow)[1:])>=dMinSDCuttOff):
				liKeepFeatures.append(iRowIndex)
		
		#Compress array
		self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]

		#Update filter state
		self._strCurrentFilterState += ":dMinSDCuttOff=" + str(dMinSDCuttOff)

		return True

        #Happy path tested 2 tests
	def funcGetWithoutOTUs(self):
		"""
		Remove features that are terminal otus. Terminal otus are identified as being an integer.
		"""

		#Get the feature names
		lsFeatures = self.funcGetFeatureNames()

		#Reduce, filter the feature names
		lsFeatures = [sFeature for sFeature in lsFeatures if not (ValidateData.funcIsValidStringInt(sFeature.split(self.funcGetFeatureDelimiter())[-1]))]

		return self.funcGetFeatureAbundanceTable(lsFeatures)

	#Happy path tested
	def funcNormalize(self):
		"""
		Convenience method which will call which ever normalization is approriate on the data.
		:return Boolean: Indicator of success (true).
		"""

		if self._fIsSummed:
			return self.funcNormalizeColumnsWithSummedClades()
		else:
			return self.funcNormalizeColumnsBySum()

	#Testing Status: Light happy path testing
	def funcNormalizeColumnsBySum(self):
		"""
		Normalize the data in a manner that is approrpiate for NOT summed data.
		Normalize the columns (samples) of the abundance table.
		Normalizes as a fraction of the total (number/(sum of all numbers in the column)).
		Will not act on summed tables.

		:return	Boolean:	Indicator of success. False indicates error.
		"""

		if self._fIsNormalized:
#			sys.stderr.write( "This table is already normalized, did not perform new normalization request.\n" )
			return False

		if self._fIsSummed:
			sys.stderr.write( "This table has clades summed, this normalization is not appropriate. Did not perform.\n" )
			return False

		#Normalize
		for columnName in self.funcGetSampleNames():
			column = self._npaFeatureAbundance[columnName]
			columnTotal = sum(column)
			if(columnTotal > 0.0):
				column = column/columnTotal
			self._npaFeatureAbundance[columnName] = column

		#Indicate normalization has occured
		self._fIsNormalized = True

		return True

	#Happy path tested
	def funcNormalizeColumnsWithSummedClades(self):
		"""
		Normalizes a summed Abundance Table.
		If this is called on a dataset which is not summed and not normalized.
		The data will be summed first and then normalized.
		If already normalized, the current normalization is kept.

		:return	Boolean:	Indicator of success. False indicates error.
		"""

		if self._fIsNormalized:
#			sys.stderr.write( "This table is already normalized, did not perform new normalization request.\n" )
			return False

		if not self._fIsSummed:
#			sys.stderr.write( "This table does not have clades summed, this normalization is not appropriate until the clades are summed. The clades are being summed now before normalization.\n" )
			self.funcSumClades()

		#Load a hash table with root data {sKey: npaAbundances}
		hashRoots = {}
		for npaRow in self._npaFeatureAbundance:

			curldAbundance = np.array(list(npaRow)[1:])
			curFeatureNameLength = len(npaRow[0].split(self._cFeatureDelimiter))
			curlRootData = hashRoots.get(npaRow[0].split(self._cFeatureDelimiter)[0])

			if not curlRootData:
				hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]] = [curFeatureNameLength, curldAbundance]
			elif curlRootData[0] > curFeatureNameLength:
				hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]] = [curFeatureNameLength, curldAbundance]

		#Normalize each feature by thier root feature
		dataMatrix = list()
		for npaRow in self._npaFeatureAbundance:

			curHashRoot = list(hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]][1])
			dataMatrix.append(tuple([npaRow[0]]+[npaRow[i+1]/curHashRoot[i] if curHashRoot[i] > 0 else 0 for i in xrange(len(curHashRoot))]))

		self._npaFeatureAbundance = np.array(dataMatrix,self._npaFeatureAbundance.dtype)

		#Indicate normalization has occured
		self._fIsNormalized = True

		return True
	
	def _funcRankAbundanceHelper( self, aaTodo, iRank, lRankAbundance ):
		"""
		Helper method for ranking abudance which are tied.

		:params aaTodo: List of tied ranks to change to a rank.
		:type:	List of Enumerates of samples.
		:params iRank: Current Rank
		:type:	Integer
		:params lRankAbundance: Sample of abundance
		:type:	List of integers
		"""

		# Subtract one from iRank (each time) to account for next loop iteration
		# Then average it with itself minus (the length of aaTodo + 1)
		dRank = ( iRank + iRank - len( aaTodo ) - 1 ) / 2.0
		for a in aaTodo:
			lRankAbundance[a[0]] = dRank

	#1 Happy path test
	def funcRankAbundance(self):
		"""
		Rank abundances of features with in a sample.

		:return	AbundanceTable:	Abundance table data ranked (Features with in samples).
							  None is returned on error.
		"""

		if self._npaFeatureAbundance == None:
			return None

		lsSampleNames = self.funcGetSampleNames()
		npRankAbundance = self.funcGetAbundanceCopy()
		liRanks = []
		#For each sample name get the ranks
		for sName in lsSampleNames:
			#Enumerate for order and sort abundances
			lfSample = list(enumerate(npRankAbundance[sName]))
			lfSample = sorted(lfSample, key = lambda a: a[1], reverse = True)

			# Accumulate indices until a new value is encountered to detect + handle ties
			aaTodo = []
			for i, a in enumerate( lfSample ):
				if ( not aaTodo ) or ( a[1] == aaTodo[-1][1] ):
					aaTodo.append( a )
				else:
			# Make multiple tied ranks = average of first and last
					self._funcRankAbundanceHelper( aaTodo, i, npRankAbundance[sName] )
					aaTodo = [a]
			self._funcRankAbundanceHelper( aaTodo, i + 1, npRankAbundance[sName] )

		return AbundanceTable(npaAbundance=npRankAbundance, dictMetadata=self.funcGetMetadataCopy(),
			strName= self.funcGetName() + "-Ranked",
			strLastMetadata=self.funcGetLastMetadataName(),
			cFileDelimiter=self.funcGetFileDelimiter(),
			cFeatureNameDelimiter=self.funcGetFeatureDelimiter())

	#Happy Path Tested
	def funcReduceFeaturesToCladeLevel(self, iCladeLevel):
		"""
		Reduce the current table to a certain clade level.

		:param	iCladeLevel:	The level of the clade to trim the features to.
		:type:	Integer	The higher the number the more clades are presevered in the consensus lineage contained in the feature name.
		:return	Boolean:	Indicator of success. False indicates error.
		"""

		if iCladeLevel < 1: return False
		if not self._npaFeatureAbundance == None:
			liFeatureKeep = []
			[liFeatureKeep.append(tplFeature[0]) if (len(tplFeature[1][0].split(self.funcGetFeatureDelimiter())) <= iCladeLevel) else 0
			 for tplFeature in enumerate(self._npaFeatureAbundance)]
			#Compress array
			self._npaFeatureAbundance = self._npaFeatureAbundance[liFeatureKeep,:]

			#Update filter state
			self._strCurrentFilterState += ":iCladeLevel=" + str(iCladeLevel)
			return True
		else:
			return False

	#Happy path tested
	def funcRemoveSamples(self,lsSampleNames):
		"""
		Removes the samples given in the list.

		:param	lsSampleNames:	A list of string names of samples to remove.
		:type:	List of strings	Unique values
		:return Boolean: Indicator of success (True = success, no error)
		"""

		#Samples to remove
		setSamples = set(lsSampleNames)

		#Get orignal sample count
		iOriginalCount  = self._iOriginalSampleCount

		#The samples to keep
		lsKeepSamples = [sSample for sSample in self.funcGetSampleNames() if not sSample in setSamples]
		#The sample to keep as boolean flags for compressing the metadata
		lfKeepSamples = [not sSample in setSamples for sSample in self.funcGetSampleNames()]
		
		#Reduce the abundance data and update
		self._npaFeatureAbundance = self._npaFeatureAbundance[[self.funcGetIDMetadataName()]+lsKeepSamples]

		#Reduce the metadata and update
		for sKey in self._dictTableMetadata:
			self._dictTableMetadata[sKey] = [value for iindex, value in enumerate(self._dictTableMetadata[sKey]) if lfKeepSamples[iindex]]

		#Update sample number count
		self._iOriginalSampleCount = len(self.funcGetSampleNames())

		return self._iOriginalSampleCount == (iOriginalCount-len(setSamples))

	#Happy path tested
	def funcRemoveSamplesByMetadata(self, sMetadata, lValuesToRemove):
		"""
		Removes samples from the abundance table based on values of a metadata.
		If a metadata has any value given the associated sample is removed.

		:param	sMetadata:	ID of the metdata to check the given values.
		:type:	String	Metadata ID
		:param	lValuesToRemove:	A list of values which if equal to a metadata entry indicate to remove the associated sample.
		:type:	List of values:	List
		:return	Boolean:	Indicator of success (True = success, no error)
		"""

		lsSampleNames = self.funcGetSampleNames()
		return self.funcRemoveSamples([lsSampleNames[iindex] for iindex, sValue in enumerate(self.funcGetMetadata(sMetadata)) if sValue in lValuesToRemove])

	#Happy path testing
	def funcSumClades(self):
		"""
		Sums abundance data by clades indicated in the feature name (as consensus lineages).

		:return	Boolean:	Indicator of success.
					False indicates an error.
		"""

		if not self.funcIsSummed():

			#Read in the data
			#Find the header column (iCol) assumed to be 1 or 2 depending on the location of "NAME"
			#Create a list (adSeq) that will eventually hold the sum of the columns of data
			astrHeaders = iCol = None
			adSeqs = np.array([0] * len(self.funcGetSampleNames()))
			pTree = CClade( )
			aastrRaw = []

			#For each row in the npaAbundance
			#Get the feature name, feature abundances, and sum up the abudance columns
			#Keep the sum for later normalization
			#Give a tree the feature name and abundance
			for dataRow in self._npaFeatureAbundance:
				
				sFeatureName = dataRow[0]
				ldAbundances = list(dataRow)[1:]

				#Add to the sum of the columns (samples)
				adSeqs = adSeqs + np.array(list(dataRow)[1:])

				#Build tree
				pTree.get( sFeatureName.split(self._cFeatureDelimiter) ).set( ldAbundances )

			#Create tree of data
			#Input missing data
			#Fill hashFeatures with the clade name (key) and a blist of values (value) of the specified level interested.
			pTree.impute( )
			hashFeatures = {}
			pTree.freeze( hashFeatures, c_iSumAllCladeLevels, c_fOutputLeavesOnly )
			setstrFeatures = hashFeatures.keys( )

			#Remove parent clades that are identical to child clades
			for strFeature, adCounts in hashFeatures.items( ):
					astrFeature = strFeature.strip( ).split( "|" )
					while len( astrFeature ) > 1:
						astrFeature = astrFeature[:-1]
						strParent = "|".join( astrFeature )
						adParent = hashFeatures.get( strParent )
						if adParent == adCounts:
							del hashFeatures[strParent]
							setstrFeatures.remove( strParent )

			#Sort features to be nice
			astrFeatures = sorted( setstrFeatures )

			#Change the hash table to an array
			dataMatrix = list()
			for sFeature in astrFeatures:
				dataMatrix.append(tuple([sFeature]+list(hashFeatures[sFeature])))
			self._npaFeatureAbundance=np.array(dataMatrix,self._npaFeatureAbundance.dtype)

			#Indicate summation has occured
			self._fIsSummed = True

		return True

	#Happy path tested
	def funcStratifyByMetadata(self, strMetadata, fWriteToFile=False):
		"""
		Stratifies the AbundanceTable by the given metadata.
		Will write each stratified abundance table to file
		if fWriteToFile is True the object will used it's internally stored name as a file to write to
		if fWriteToFile is a string then it should be a directory and end with "." This will rebase the file
		and store it in a different directory but with an otherwise unchanged name.
		Note: If the metadata used for stratification has NAs, they will be segregated to thier own table and returned.

		:param	strMetadata:	Metadata ID to stratify data with.
		:type:	String	ID for a metadata.
		:param	fWriteToFile:	Indicator to write to file.
		:type:	Boolean	True indicates to write to file.
		:return	List:	List of AbundanceTables which are deep copies of the original.
						Empty list on error.
		"""

		if self._npaFeatureAbundance is None or self._dictTableMetadata is None:
			return []

		#Get unique metadata values to stratify by
		lsMetadata = self._dictTableMetadata.get(strMetadata,[])
		setValues = set(lsMetadata)
		#If there is only one metadata value then no need to stratify so return the original in the list (and write if needed)
		if len(setValues) == 0:
		  return []

		retlAbundanceTables = []
		dictAbundanceBlocks = dict()
		#Given here there are multiple metadata values, continue to stratify
		lsNames = self.funcGetSampleNames()
		#Get index of values to break up
		for value in setValues:
			lfDataIndex = [sData==value for sData in lsMetadata]
			#Get abundance data for the metadata value
			#The true is added to keep the first column which should be the feature id
			npaStratfiedAbundance = self._npaFeatureAbundance[[self.funcGetIDMetadataName()]+list(np.compress(lfDataIndex,lsNames))]

			#Get metadata for the metadata value
			dictStratifiedMetadata = dict()
			for metadataType in self._dictTableMetadata:
				dictValues = self.funcGetMetadata(metadataType)
				dictStratifiedMetadata[metadataType] = np.compress(lfDataIndex,dictValues).tolist()

			#Make abundance table
			#Add abundance table to the list
			lsNamePieces = os.path.splitext(self._strOriginalName)
			objStratifiedAbundanceTable = AbundanceTable(npaAbundance=npaStratfiedAbundance, dictMetadata=dictStratifiedMetadata,
				strName=lsNamePieces[0] + "-StratBy-" + value+lsNamePieces[1],
				strLastMetadata=self.funcGetLastMetadataName(),
				cFeatureNameDelimiter=self._cFeatureDelimiter, cFileDelimiter = self._cDelimiter)
			if fWriteToFile:
				objStratifiedAbundanceTable.funcWriteToFile(lsNamePieces[0] + "-StratBy-" + value+lsNamePieces[1])
			#Append abundance table to returning list
			retlAbundanceTables.append(objStratifiedAbundanceTable)

		return retlAbundanceTables

	#Happy Path Tested
	def funcTranslateIntoMetadata(self, lsValues, sMetadataFrom, sMetadataTo, fFromPrimaryIds=True):
		"""
		Takes the given data values in one metadata and translates it to values in another
		metadata of the sample samples holding the values of the first metadata
		FPrimaryIds, if true the sMetadataFrom are checked for unique values,
		If FPrimaryIds is not true, duplicate values can stop the preservation of order
		Or may cause duplication in the "to" group. This is not advised.
		if the sMetadataFrom has any duplicates the function fails and return false.

		:param	lsValues:	Values to translate.
		:type:	List	List of values.
		:param	sMetadataFrom:	The metadata the lsValues come from.
		:type:	String	ID for the metadata.
		:param	sMetadataTo:	The metadata the lsValues will be translated into keeping the samples the same.
		:type:	String	ID for the metadata.
		:param	fFromPrimaryIds:	The metadata that are in the from metadata list must be unique in each sample.
		:type:	Boolean	True indicates the metadata list should be unique in each sample. Otherwise a false will return.
		:return List:	List of new values or False on error.
		"""

		#Get metadata
		lFromMetadata = self.funcGetMetadata(sMetadataFrom)
		if not lFromMetadata:
				sys.stderr.write( "Abundancetable::funcTranlateIntoMetadata. Did not receive lFromMetadata.\n" )
				return False

		lToMetadata = self.funcGetMetadata(sMetadataTo)
		if not lToMetadata:
				sys.stderr.write( "Abundancetable::funcTranlateIntoMetadata. Did not receive lToMetadata.\n" )
				return False

		#Check to see if the values are unique if indicated to do so
		if fFromPrimaryIds:
			if not len(lFromMetadata) == len(set(lFromMetadata)):
				sys.stderr.write( "Abundancetable::funcTranlateIntoMetadata. sMetadataFrom did not have unique values.\n" )
				return False

		#Translate over
		if lFromMetadata and lToMetadata:
			return [lToMetadata[iIndex] for iIndex in [lFromMetadata.index(value) for value in lsValues]]

		return False

	#Happy path tested
	def funcToArray(self):
		"""
		Returns a numpy array of the current Abundance Table.
		Removes the first ID head column and the numpy array is
		Made of lists, not tuples.

		:return Numpy Array:	np.array([[float,float,...],[float,float,...],[float,float,...]])
								None is returned on error.
		"""

		if not self._npaFeatureAbundance == None:
			return np.array([list(tplRow)[1:] for tplRow in self._npaFeatureAbundance],'float')
		return None

	#Happy Path tested
	def funcWriteToFile(self, xOutputFile, cDelimiter=None):
		"""
		Writes the AbundanceTable to a file strOutputFile.
		Will rewrite over a file as needed.
		Will use the cDelimiter to delimit columns if provided.

		:param	xOutputFile:	File stream or File path to write the file to.
		:type:	String	File Path
		:param	cDelimiter:	Delimiter for the output file.
		:type:	Character	If cDlimiter is not specified, the internally stored file delimiter is used.
		"""

		if not xOutputFile:
			return
		#Check delimiter argument
		if not cDelimiter:
			cDelimiter = self._cDelimiter

		f = csv.writer(open( xOutputFile, "w" ) if isinstance(xOutputFile, str) else xOutputFile, csv.excel_tab, delimiter=cDelimiter)
		
		#Write Ids
		f.writerows([[self.funcGetIDMetadataName()]+list(self.funcGetSampleNames())])
		#Write metadata
		lsKeys = list(set(self._dictTableMetadata.keys())-set([self.funcGetIDMetadataName(),self.funcGetLastMetadataName()]))
		f.writerows([[sMetaKey]+self.funcGetMetadata(sMetaKey) for sMetaKey in lsKeys+[self.funcGetLastMetadataName()]])
		#Write abundance
		lsOutput = list()
		curAbundance = self._npaFeatureAbundance.tolist()
		for curAbundanceRow in curAbundance:
			f.writerows([[str(curAbundanceElement) for curAbundanceElement in curAbundanceRow]])

	#Testing Status: 1 Happy path test
	@staticmethod
	def funcPairTables(strFileOne, strFileTwo, strIdentifier, cDelimiter, strOutFileOne, strOutFileTwo, lsIgnoreValues=None):
		"""
		This method will read in two files and abridge both files (saved as new files)
		to just the samples in common between the two files given a common identifier.
		***If the identifier is not unique in each data set, the first sample with the pairing id is taken so make sure the ID is unique.
		Expects the files to have the sample delimiters.

		:param	strFileOne:	Path to file one to be paired.
		:type:	String	File path.
		:param	strFileTwo:	Path to file two to be paired.
		:type:	String	File path.
		:param	strIdentifier:	Metadata ID that is used for pairing.
		:type:	String	Metadata ID.
		:param	cDelimiter:	Character delimiter to read the files.
		:type:	Character	Delimiter.
		:param	strOutFileOne:	The output file for the paired version of the first file.
		:type:	String	File path.
		:param	strOutFileTwo:	The output file for the paired version of the second file.
		:type:	String	File path.
		:param	lsIgnoreValues:	These values are ignored even if common IDs between the two files.
		:type:	List	List of strings.
		:return	Boolean:	Indicator of no errors.
							  False indicates errors.
		"""

		#Validate parameters
		if(not ValidateData.funcIsValidFileName(strFileOne)):
			sys.stderr.write( "AbundanceTable:checkRawDataFile::Error, file not valid. File:" + strFileOne + "\n" )
			return False
		#Validate parameters
		if(not ValidateData.funcIsValidFileName(strFileTwo)):
			sys.stderr.write( "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+ strFileTwo + "\n" )
			return False

		#Make file one
		#Read in file
		istm = csv.reader(open(strFileOne,'r'), csv.excel_tab, delimiter=cDelimiter)
		lsContentsOne = [lsRow for lsRow in istm]

		#Get the file identifier for file one
		fileOneIdentifier = None
		for sLine in lsContentsOne:
			if sLine[0] == strIdentifier:
				fileOneIdentifier = sLine
				break

		#Make file two
		#Read in file
		istm = csv.reader(open(strFileTwo,'r'), csv.excel_tab, delimiter=cDelimiter)
		lsContentsTwo = [lsRow for lsRow in istm]

		#Get the file identifier for file two
		fileTwoIdentifier = None
		for sLine in lsContentsTwo:
			if sLine[0] == strIdentifier:
				fileTwoIdentifier = sLine
				break

		#Get what is in common between the identifiers
		#And find which columns to keep in the tables based on the common elements
		setsCommonIdentifiers = set(fileOneIdentifier) & set(fileTwoIdentifier)
		if lsIgnoreValues:
			setsCommonIdentifiers = setsCommonIdentifiers - set(lsIgnoreValues)

		#Get positions of common identifiers in each data set, if the identifier is not unique in a date set just take the first index
		lfFileOneIDIndexes = [fileOneIdentifier.index(sCommonID) for sCommonID in setsCommonIdentifiers]
		lfFileTwoIDIndexes = [fileTwoIdentifier.index(sCommonID) for sCommonID in setsCommonIdentifiers]

		#Convert index list to list of boolean
		lfFileOneElements = [iIndex in lfFileOneIDIndexes for iIndex, sIdentifier in enumerate(fileOneIdentifier)]
		lfFileTwoElements = [iIndex in lfFileTwoIDIndexes for iIndex, sIdentifier in enumerate(fileTwoIdentifier)]

		#Write out file one
		ostm = csv.writer(open(strOutFileOne,'w'), csv.excel_tab, delimiter=cDelimiter)
		(ostm.writerows([np.compress(lfFileOneElements,sLine) for sLine in lsContentsOne]))

		#Write out file two
		ostm = csv.writer(open(strOutFileTwo,'w'), csv.excel_tab, delimiter=cDelimiter)
		(ostm.writerows([np.compress(lfFileTwoElements,sLine) for sLine in lsContentsTwo]))

		return True

	#Testing Status: Light happy path testing
	@staticmethod
	def funcStratifyAbundanceTableByMetadata(strInputFile = None, strDirectory = "", cDelimiter = ConstantsBreadCrumbs.c_cTab, iStratifyByRow = 1, llsGroupings = []):
		"""
		Splits an abundance table into multiple abundance tables stratified by the metadata

		:param	strInputFile:	String file path to read in and stratify.
		:type:	String	File path.
		:param	strDirectory:	Output directory to write stratified files.
		:type:	String	Output directory path.
		:param	cDelimiter:	The delimiter used in the adundance file.
		:type:	Character	Delimiter.
		:param	iStratifyByRow:	The row which contains the metadata to use in stratification.
		:type:	Integer	Positive integer index.
		:param	llsGroupings:	A list of string lists where each string list holds values that are equal and should be grouped together.
								So for example, if you wanted to group metadata "1", "2", and "3" seperately but "4" and "5" together you would
								Give the following [["4","5"]].
								If you know what "1" and "3" also together you would give [["1","3"],["4","5"]]
		:type	List	List of list of strings
		:return	Boolean:	Indicator of NO error.
							False indicates an error.
		"""

		#Validate parameters
		if(not ValidateData.funcIsValidFileName(strInputFile)):
			sys.stderr.write( "AbundanceTable:stratifyAbundanceTableByMetadata::Error, file not valid. File:" + strInputFile + "\n" )
			return False
		if(not ValidateData.funcIsValidStringType(cDelimiter)):
			sys.stderr.write( "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Delimiter is not a valid string/char type. Delimiter =" + cDelimiter + "\n" )
			return False
		if(not ValidateData.funcIsValidPositiveInteger(iStratifyByRow, tempZero = True) and (not ValidateData.funcIsValidString(iStratifyByRow))):
			sys.stderr.write( "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Stratify by row is not a positive integer or string keyword. Row =" +
				str(iStratifyByRow) + ".\n" )
			return False

		#Get the base of the file path
		#This is dependent on the given output directory and the prefix of the file name of the input file
		#If no output file is given then the input file directory is used.
		baseFilePath = strDirectory
		lsFilePiecesExt = os.path.splitext(strInputFile)
		if baseFilePath:
			baseFilePath = baseFilePath + os.path.splitext(os.path.split(strInputFile)[1])[0]
		else:
			baseFilePath = lsFilePiecesExt[0]

		#Read in file
		istm = csv.reader(open(strInputFile,'r'), csv.excel_tab, delimiter=cDelimiter)
		sFileContents = [lsRow for lsRow in istm]

		#Collect metadata
		metadataInformation = dict()

		#If the tempStratifyRow is by key word than find the index
		if ValidateData.funcIsValidString(iStratifyByRow):
			for iLineIndex, strLine in enumerate(sFileContents):
				if strLine[0].strip("\"") == iStratifyByRow:
					iStratifyByRow = iLineIndex
					break

		#Stratify by metadata row
		#Split metadata row into metadata entries
		#And put in a dictionary containing {"variable":[1,2,3,4 column index]}
		stratifyByRow = sFileContents[iStratifyByRow]
		for metaDataIndex in xrange(1,len(stratifyByRow)):
			metadata = stratifyByRow[metaDataIndex]
			#Put all wierd categories, none, whitespace, blank space metadata cases into one bin
			if not metadata or metadata in string.whitespace:
				metadata = "Blank"
			#Remove any extraneous formatting
			metadata = metadata.strip(string.whitespace)
			#Store processed metadata with column occurence in dictionary
			if(not metadata in metadataInformation):
				metadataInformation[metadata] = []
			metadataInformation[metadata].append(metaDataIndex)

		#For each of the groupings
		#Use the first value as the primary value which the rest of the values in the list are placed into
		#Go through the dict holding the indices and extend the list for the primary value with the secondary values
		#Then set the secondary value list to empty so that it will be ignored.
		if llsGroupings:
			for lSKeyGroups in llsGroupings:
				if len(lSKeyGroups) > 1:
					for sGroup in lSKeyGroups[1:]:
						if sGroup in metadataInformation:
							metadataInformation[lSKeyGroups[0]].extend(metadataInformation[sGroup])
							metadataInformation[sGroup] = []

		#Stratify data
		stratifiedAbundanceTables = dict()
		for tableRow in sFileContents:
			if(len(tableRow)> 1):
				for metadata in metadataInformation:
					#[0] includes the taxa line
					columns = metadataInformation[metadata]
					if columns:
						columns = [0] + columns
						lineList = list()
						for column in columns:
							lineList.append(tableRow[column])
						stratifiedAbundanceTables.setdefault(metadata,[]).append(lineList)

		#Write to file
		lsFilesWritten = []
		for metadata in stratifiedAbundanceTables:
			sOutputFile = baseFilePath+"-by-"+metadata.strip("\"")+lsFilePiecesExt[1]
			f = csv.writer(open(sOutputFile,'w'), csv.excel_tab, delimiter = cDelimiter )
			f.writerows(stratifiedAbundanceTables[metadata])
			lsFilesWritten.append(sOutputFile)

		return lsFilesWritten
