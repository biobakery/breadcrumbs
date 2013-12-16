# BreadCrumbs #

BreadCrumbs is an unofficial collection of scripts and code intended to consolidate functions for tool development and contain scripts for command line access to commonly used functions. Breadcrumbs tends to include functionality associated with metagenomics analysis but you never know what you will find!


## Dependencies: ##

NOTE: Dependencies will be installed, if needed, with the setup.py install function below.
Manual installs of the following dependencies should not be required.

1. Cogent 		https://pypi.python.org/pypi/cogent
2. MatplotLib 		http://matplotlib.org/downloads.html
3. Mercurial 		http://mercurial.selenic.com/ (optional for downloading)
4. Numpy 		http://www.numpy.org/
5. Python 2.x 		http://www.python.org/download/
6. SciPy 		http://www.scipy.org/install.html
7. biom support 	http://biom-format.org/
8. blist		https://pypi.python.org/pypi/blist/
9. R			http://www.r-project.org/
10. R libraries including vegan, optparse

## How to download ##

To download BreadCrumbs from BitBucket use the command:

> hg clone https://bitbucket.org/timothyltickle/breadcrumbs

To update BreadCrumbs, in the BreadCrumbs directory use the 2 commands sequentially:

> hg pull  
> hg update


## Install: ##

To install please use the setup.py script in a terminal.

> python setup.py install


## Note to Mac users: ##

Before using the setup.py script you will need to install:
GCC compiler		Go to https://developer.apple.com/downloads/index.action
			Register for an Apple account you do not have one (it is free)
			In the Developer Tools category download "Command Line Tools for Xcode"
			Download and install.
Fortran Compiler	See Fortran section of http://www.scipy.org/scipylib/building/macosx.html

You may also need to run the command: (errors referencing the egg).
> pip install setuptools --upgrade


## Scripts: ##

Scripts are included to expose core functionality through the command line. Currently these scripts center on manipulating and visualizing abundance tables.  
A quick description of the scripts include:

* *Hclust.py* Flexible script to create a visualization of hierarchical clustering of abundance tables (or other matrices).

* *scriptBiplotTSV.R* Allows one to plot a tsv file as a biplot using nonmetric multidimensional scaling.

* *scriptPlotFeature.py* Allows one to plot a histogram, boxplot, or scatter plot of a bug or metadata in an abundance table. Will work on any row in a matrix.

* *scriptManipulateTable.py* Allows one to perform common functions on an abundance table including, summing, normalizing, filtering, stratifying tables.

* *scriptPcoa.py* Allows one to plot a principle covariance analysis (PCoA) plot of an abundance table.

* *scriptConvertBetweenBIOMAndPCL.py* Allows one to convert between BIOM and PCL file formats.


## Programming Classes: ##

Brief descriptions of classes are as follows. More detailed descriptions are given in the classes themselves.

* *AbundanceTable* Data structure to contain and perform operations on an abundance table.

* *BoxPlot* Wrapper to plot box plots.

* *CClade* Helper object used in hierarchical summing and normalization

* *Cladogram* Object that manipulated an early dendrogram visualization. Deprecated, should use the GraPhlan visualization tool on bitbucket instead.

* *CommandLine* Collection of code to work with command line. Deprecated. Should use sfle calls.

* *ConstantsBreadCrumbs* Contains generic constants.

* *ConstantsFiguresBreadCrumbs* Contains constants associated with formatting figures.

* *KMedoids* Code from MLPY which performs KMedoids sample selection.

* *MLPYDistanceAdaptor* Used to allow custom distance matrices to be used by KMedoids.

* *Metric* Difference functions associated with distance and diversity metrics.

* *PCoA* Functionality surrounding the plotting of a PCoA

* *PlotMatrix* Allows on to plot a matrix of numbers.

* *SVM* Support Vector Machine associated scripts.

* *Utility* Generic functions

* *UtilityMath* Generic math related functions

* *ValidateData* Collection of functions to validate data types when needed.


## Demo input files: ##

* *fastunifrac_Ley_et_al_NRM_2_sample_id_map.txt* Example Unifrac Id mapping file (source http://bmf2.colorado.edu/fastunifrac/tutorial.psp)

* *GreenGenesCore-May09.ref.tre* Example Greengenes core set reference for Unifrac demo (source http://bmf2.colorado.edu/fastunifrac/tutorial.psp)

* *Test.pcl* Example file / Test PCL file to run scripts on.

* *Test.biom* Example file / Test BIOM file to run scripts on.

* *Test_no_metadata.pcl* Example file / Test PCL file to run scripts on which does not have metadata.

* *Test_no_metadata.biom* Example file / Test BIOM file to run scripts on which does not have metadata.

* *Test-biplot.tsv* Example file / Test file for the scriptBiplotTSV.R


## Contributing Authors: ##
Timothy Tickle, George Weingart, Nicola Segata, Randall Schwager, Curtis Huttenhower


## Contact: ##
Please feel free to contact ttickle@hsph.harvard.edu with questions.
