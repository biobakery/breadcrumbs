# BreadCrumbs #

BreadCrumbs is an unofficial collection of scripts and code intended to consolidate functions for tool development and contain scripts for command line access to commonly used function. Breadcrumbs tends to include functionality associated with metagenomics analysis but you never know what you will find!


## Dependencies: ##

Cogent (https://pypi.python.org/pypi/cogent)
MatplotLib (http://matplotlib.org/downloads.html)
Mercurial (optional for downloading) (http://mercurial.selenic.com/)
Numpy (http://www.numpy.org/)
Python 2.x (http://www.python.org/download/)
SciPy (http://www.scipy.org/install.html)
biome support (http://biom-format.org/)


## How to download ##

To download BreadCrumbs from BitBucket use the command:

> hg clone https://bitbucket.org/timothyltickle/breadcrumbs

To update BreadCrumbs, in the BreadCrumbs directory use the 2 commands sequentially:

> hg pull  
> hg update


## Scripts: ##

Scripts are included to expose core functionality to the command line. Currently these scripts center on manipulating and visualizing abundance tables.  
A quick description of the scripts include:

* *Hclust.py* Flexible script to create a visualization of hierarchical clustering of abundance tables (or other matrices).

* *scriptBiplotPCL.R* Allows one to plot a pcl as a biplot using nonmetric multidimensional scaling.

* *scriptPlotFeature.py* Allows one to plot a histogram, boxplot, or scatter plot of a bug or metadata in an abundance table. Will work on any row in a matrix.

* *scriptManipulateTable.py* Allows one to do common functions on an abundance table including, summing, normalizing, filtering, stratifying tables.

* *scriptPcoa.py* Allows one to plot a principle covariance analysis (PCoA) plot of an abundance table.

* *scriptConvertBetweenBiomeAndPCL.py* Allows one to convert between BIOME and PCL file formats.


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

* *Test.biome* Example file / Test BIOME file to run scripts on.

* *Test_no_metadata.pcl* Example file / Test PCL file to run scripts on which does not have metadata.

* *Test_no_metadata.biome* Example file / Test BIOME file to run scripts on which does not have metadata.

* *Test-biplot.tsv* Example file / Test file for the scriptBiplotTSV.R


## Contributing Authors: ##
Timothy Tickle, George Weingart, Nicola Segata, Curtis Huttenhower


## Contact: ##
Please feel free to contact ttickle@hsph.harvard.edu with questions.
