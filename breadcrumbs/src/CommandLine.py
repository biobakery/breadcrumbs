"""
Author: Timothy Tickle
Description: Manages calling commandline from within code.
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

#Import libaries
from subprocess import call, Popen, PIPE
from ValidateData import ValidateData
import traceback

class CommandLine():
    """
    Manages calling commandline from within code.
    """

    ##
    #Contructor
    def __init__(self): pass

    def runCommandLine(self,tempCommand = None):
        """
        Sends a command to command line interface.

        :param tempCommand: Must be an list of command key word and string arguments, no whitespaces
        :type: List of strings
        :return: boolean indicator of success (True = Success)
        """

        #Makes sure the the input data is a list of strings
        if(not ValidateData.funcIsValidStringList(tempCommand)):
            print "Error:: tempCommand must be an array of strings. Received="+str(tempCommand)+"."
            return False

        #Run command
        try:
            returnCode = call(tempCommand)
            print "Return="+str(returnCode)
            if returnCode > 0:
                print "Error:: Error during command call. Script stopped."
                print "Error:: Error Code "+str(returnCode)+"."
                print "Error:: Command ="+str(tempCommand)+"."
                return False
        except (OSError,TypeError), e: 
                print "Error:: Error during command call. Script stopped."
                print "Error:: Command ="+str(tempCommand)+"."
                print "Error:: OS error: "+str(traceback.format_exc(e))+"."
                return False
        return True

    def runPipedCommandLine(self,tempCommand = None):
        """
        Sends a command to command line interface.
        Create new array of string elements instead of white spacing
        Put file names in escaped quotation marks.
        This uses shell == true so make sure the commandline is not malicious
        This should wait for process completion

        :param tempCommand: Must be an list of command key word and string arguments, no whitespaces.
        :type: List of strings
        :return: Boolean (False = Failure or the return code from the subprocess)
        """

        #Makes sure the the input data is a list of strings
        if(not ValidateData.funcIsValidStringList(tempCommand)):
            print "Error:: tempCommand must be an array of strings. Received="+str(tempCommand)+"."
            return False

        #Run command
        tempCommand = " ".join(tempCommand)
        try:
            returnCode = Popen(tempCommand, shell = True, stdout = PIPE).communicate()
            return returnCode
        except (OSError,TypeError), e: 
                print "Error:: Error during command call. Script stopped."
                print "Error:: Command ="+str(tempCommand)+"."
                print "Error:: OS error: "+str(traceback.format_exc(e))+"."
                return False

    def runBatchCommandline(self,tempArrayOfCommands = None):
        """
        Sends a an array of commands to the commandline.

        :param tempArrayOfCommands: Must be an list of commands, parsing and removing whitespace is handled internally.
         Do not send mkdir and rm commands, use the appropriate os.* method call
        :type: List of strings
        :return: boolean indicator of success (True = Success)
        """

        #Holds commands
        parsedCommmands = []

        #Indicates if success or error occured
        success = True

        #Makes sure the the input data is list of strings
        if(not ValidateData.funcIsValidStringList(tempArrayOfCommands)):
            print "Error:: tempCommand must be an array of strings. Received="+str(tempArrayOfCommands)+"."
            return False

        #Parse commands into an array and call
        #On an error break and return False
        for command in tempArrayOfCommands:
            commandElements = command.split(" ")
            if(not self.runCommandLine(commandElements)):
                success = False
                break
        return success
                
