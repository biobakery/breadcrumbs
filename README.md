# BreadCrumbs #

BreadCrumbs is an unofficial collection of scripts and code intended to consolidate functions for tool development and contain scripts for command line access to core function. Breadcrumbs tends to include functionality associated with metagenomics analysis but you never know what you will find!


## Dependencies: ##

Cogent
MatplotLib
Mercurial (optional for downloading)
Numpy
Python 2.x
SciPy


## How to download ##

To download BreadCrumbs from BitBucket use the command:

> hg clone https://bitbucket.org/timothyltickle/breadcrumbs

To update BreadCrumbs, in the BreadCrumbs directory use the 2 commands sequentially:

> hg pull
> hg update


## Scripts: ##

Scripts are included to expose core functionality to the command line. Currently these scripts center on Manipulating and visualizing abundance tables.
A quick description of the scripts include:

* Hclust.py
Flexible script to create hierarchical clustering of abundance tables (or other matrices).

* scriptPlotFeature.py
Allows one to plot a histogram, boxplot, or scatter plot of a bug or metadata in an abundance table. Will work on any row in a matrix.

* scriptManipulateTable.py
Allows one to do common functions on an abundance table including, summing, normalizing, filtering, stratifying tables.

* scriptPcoa.py
Allows one to plot a principle covariance analysis (PCoA) plot of an abundance table.


## Programming Classes: ##

Brief descriptions of classes are as follows. More detailed descriptions are given in the classes themselves.

* AbundanceTable
Data structure to contain and perform operations on an abundance table.

* BoxPlot
Wrapper to plot box plots.

* CClade
Helper object used in hierarchical summing and normalization

* Cladogram
Object that manipulated an early dendrogram visualization. Deprecated, should use the GraPhlan visualization tool on bitbucket instead.

* CommandLine
Collection of code to work with command line. Deprecated. Should use sfle calls.

* ConstantsBreadCrumbs
Contains generic constants.

* ConstantsFiguresBreadCrumbs
Contains constants associated with formatting figures.

* KMedoids
Code from MLPY which performs KMedoids sample selection.

* MLPYDistanceAdaptor
Used to allow custom distance matrices to be used by KMedoids.

* Metric
Difference functions associated with distance and diversity metrics.

* PCoA
Functionality surrounding the plotting of a PCoA

* PlotMatrix
Allows on to plot a matrix of numbers.

* SVM
Support Vector Machine associated scripts.

* Test.pcl
Test file to run script on.

* Utility
Generic functions

* UtilityMath
Generic math related functions

* ValidateData
Collection of functions to validate data types when needed.


## Contributing Authors: ##
Curtis Huttenhower, Nicola Segata, Timothy Tickle


## Contact: ##
Please feel free to contact ttickle@hsph.harvard.edu with questions.