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

MIN_SETUPTOOLS = '0.6.28'

import sys
import subprocess

import setuptools
from distutils.version import StrictVersion

arguments = sys.argv[1:]

if StrictVersion(setuptools.__version__) < StrictVersion(MIN_SETUPTOOLS):
    print "Your version of distribute (%s) it too low. (must be >= %s)" %(
        setuptools.__version__, MIN_SETUPTOOLS)
    print "Now upgrading distribute..."

    subprocess.check_call('pip install --upgrade distribute'.split())

    print "Complete. New version: %s" %(setuptools.__version__)

subprocess.check_call('pip install numpy==1.7.1'.split())
subprocess.check_call('pip install scipy==0.12.0'.split())
subprocess.check_call('pip install https://github.com/bipy/pyqi/archive/v0.2.0.tar.gz#egg=pyqi-0.2.0'.split())
cmd = 'python actually_setup.py'.split()
cmd.extend(arguments)
subprocess.check_call(cmd)
