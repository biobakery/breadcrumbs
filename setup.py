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

import subprocess

from setuptools import setup, find_packages

subprocess.check_call('pip install numpy==1.7.1'.split())
subprocess.check_call('pip install scipy==0.12.0'.split())
subprocess.check_call('pip install https://github.com/bipy/pyqi/archive/v0.2.0.tar.gz#egg=pyqi-0.2.0'.split())

setup(
    name='breadcrumbs',
    author="Timothy Tickle",
    author_email="ttickle@sph.harvard.edu",
    license="MIT",
    version='0.8.0',
    description='Assorted metagenomics tools',
    packages=find_packages(exclude=['ez_setup', 'tests', 'tests.*']),
    zip_safe=False,
    install_requires=[
        'Cython==0.19.1',
        'argparse',
        'cogent==1.5.3',
        'matplotlib==1.3.0',
        'numpy==1.7.1',
        'scipy==0.12.0',

        'pyqi==0.2.0',
        'biom-format==1.2.0',
    ],
    dependency_links=[
        'https://github.com/biom-format/biom-format/archive/v1.2.0.tar.gz#egg=biom-format-1.2.0',
    ],
    classifiers=[
        "Development Status :: 4 - Beta"
    ],
)
