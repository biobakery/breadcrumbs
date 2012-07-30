"""
Author: Timothy Tickle
Description: Utility class for generic functions.
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

class Utility():
    """
    Class to perform misc methods.
    """

    #Tested 6
    @staticmethod
    def getIndices(aList, dataElement):
        """
        Returns the index or indicies of the element in the list.

        :param aList: List ot search for element.
        :type List.
        :param dataElement: Element for which to search.
        :type Object of the same type as is found in the list.
        :return: List of indicies indicating where the element occurs in the list. Returns [] when the element is not in the list.
        """

        aretIndices = []
        for dataIndex in xrange(0,len(aList)):
            if(aList[dataIndex] == dataElement):
                aretIndices.append(dataIndex)
        return aretIndices

    #Tested 6
    @staticmethod
    def reduceList(aList, dataIndicies):
        """
        Reduces a list to just the data indicies given.

        :param aList: List to reduce.
        :type List
        :param dataIndicies: list of indicies to keep.
        :type List of integers
        :return: Reduced list.  Returns [] when the and empty index list is given.
        """
        return [aList[dataIndicies[dataIndex]] for dataIndex in xrange(0,len(dataIndicies))]

    #Tested 8
    @staticmethod
    def RGBToHex(adColor):
        """
        Change a RGB float to hex.

        :param adColor: A list of 3 elements which are floats between 0.0 and 1.0
        :type A list of floats
        :return: A string (HEX formatted) representation of the RGB color
        """

        charR = (hex(int(adColor[0]*255)))[2:]
        if(str(charR) == "0"):
            charR = "00"
        charG = (hex(int(adColor[1]*255)))[2:]
        if(str(charG) == "0"):
            charG = "00"
        charB = (hex(int(adColor[2]*255)))[2:]
        if(str(charB) == "0"):
            charB = "00"
        return "".join(["#",charR, charG, charB])
