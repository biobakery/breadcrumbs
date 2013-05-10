# BreadCrumbs Tutorial #

This is a brief tutorial to get you acquainted with the scripts provided in breadcrumbs. This is broke up by script and task. Examples are given using the Test.pcl file which is included in the BreadCrumbs package. Each of these commands should work from the command line in the breadcrumbs directory.

Enjoy and happy researching!

## Contents: ##
1. scriptPCoA  
2. scriptManipulateFeature.py  
I. Manipulating the measurements  
II. Filtering  
III. Filtering with knowledge of feature hierarchical relationship
IV. Manipulate samples by metadata
V. Manipulate the feature names
VI. Dimensionality Reduction
3. scriptPlotFeature.py  

## scriptPCoA.py ##
This script allows one to plot a PCoA of an abundance table. In the plot each sample is one marker. The marker shape and color is determined by a metadata (of your choice). The distances between each sample is determined by a specific beta-diversity distance metric. By default Bray-curtis distance is used. This can be changed as needed. You will notice for every call you must give it the sample id (-i) and the last metadata which should be the row before your first data (-l). This helps the scripts understand what is a data measurement and what is a metadata.

*Note this tutorial assumes you are in the breadcrumbs directory*

A. How do I make a PCoA of an abundance table, painting (coloring) it by a specific metadata?

> python scripts/scriptPcoa.py -i TID -l STSite -p STSite demo_input/Test.pcl

B. How do I make a series of PCoAs of an abundance table, one PCoA for every metadata?

If nothing is specified with -p then all metadata are painted. Note there are a max of 9 shapes to use, a metadata will be skipped if it has more than 9 levels (specific values which can be used many times). Don't worry, the script will let you know if this happens and will just skip to the next metadata.

> python scripts/scriptPcoa.py -i TID -l STSite demo_input/Test.pcl

C. How do I use a different beta-diversity distance metric instead of Bray-curtis distance?
Note the current two options are SPEARMAN and BRAY_CURTIS (default) more can be added as needed.

> python scripts/scriptPcoa.py -i TID -l STSite -m SPEARMAN demo_input/Test.pcl

D. How do I get the coordinates of the points in the PCoA plot? Use -C and give a file path to which to write.

> python scripts/scriptPcoa.py -i TID -l STSite -C coordinates.txt demo_input/Test.pcl

E. How do I get the distance matrix represented by the PCoA plot? Use -D and give a file path to which to write.

> python scripts/scriptPcoa.py -i TID -l STSite -D distances.txt demo_input/Test.pcl

## scriptManipulateFeature.py ##
Abundance tables can be difficult to manipulate. This script captures frequent tasks that may be important to manipulating an abundance table including normalization, summing, filtering, stratifying the tables into subsets (for instance breaking up a large HMP table into tables, one for each body site), and other functionality. You will notice for every call you must give it the sample id (-i) and the last metadata which should be the row before your first data (-l). This helps the scripts understand what is a data measurement and what is a metadata.

_Remember you can do multiple tasks or use multiple arguments at the same time._
_Here is an example of summing, normalizing, adding on clade prefixes, and stratifies the tables based on the STSite metadata_

> python scripts/scriptManipulateTable.py -i TID -l STSite -s -n -x -y STSite demo_input/Test.pcl

Please look at the detailed description of normalization and summation for a clear understanding of how the data is being manipulated.

*Manipulating the measurements*  
A. How do I sum a table based on clade names?

> python scripts/scriptManipulateTable.py -i TID -l STSite -s demo_input/Test.pcl

B. How do I normalize a table?

> python scripts/scriptManipulateTable.py -i TID -l STSite -n demo_input/Test.pcl

*Filtering*  
C. How do I filter a normalized table by percentage?

This filters out bugs that are not in the top 0.95 percentage of at least 0.05 percent of the samples (a good default).

> python scripts/scriptManipulateTable.py -i TID -l STSite -P 0.95,0.05 demo_input/Test.pcl

D. How do I filter a normalized table by a minimum abundance?  

This filters out bugs that do not have at least a certain number of bugs in a certain number of samples. Here we show
the call to filter out all bugs which do not have at least 3 samples with at least 0.0001 abundance (a good initial default).

> python scripts/scriptManipulateTable.py -i TID -l STSite -A 0.0001,3 demo_input/Test.pcl

E. How do I filter a count table by count occurrence?  

This removes samples that do not have at least 5 counts in at least 3 samples (an initial default to use could be 2,2).

> python scripts/scriptManipulateTable.py -i TID -l STSite -O 5,3 demo_input/Test.pcl

F. How do I filter a table by standard deviation?  

> python scripts/scriptManipulateTable.py -i TID -l STSite -D 1 demo_input/Test.pcl

*Filtering with knowledge of feature hierarchical relationship.*  
F. How do I make the table have only terminal nodes?

> python scripts/scriptManipulateTable.py -i TID -l STSite -t demo_input/Test.pcl

G. How do I remove all the OTUs from a table?

> python scripts/scriptManipulateTable.py -i TID -l STSite -u demo_input/Test.pcl

H. How do I reduce all bugs more specific than a certain clade? Aka, how do I reset a table to be only a clade (genus) or higher?  

This reduces the bugs to bugs with 3 levels of hierarchy or less (class on a standard biological taxonomy).

> python scripts/scriptManipulateTable.py -i TID -l STSite -c 3 demo_input/Test.pcl

You may want to hierarchically sum all of you bugs before reducing the table to a certain level, just in case you are missing some.

> python scripts/scriptManipulateTable.py -i TID -l STSite -s -c 3 demo_input/Test.pcl

Detail. OTUs or taxonomic clades are terminal nodes of a dendrogram representing the full taxonomy or phylogeny of a study. Biology may happen at these terminal clades or at higher level clades. Hierarchical summation uses the name of the bug (containing the consensus lineage) to add bugs together at different levels of their ancestral state and represent additional higher level clades or bigger groupings of bugs.

More plainly, imagine if we have 2 bugs in a sample with 5 and 10 counts. These two bugs differ as species but share the rest of their ancestry. In this case, an additional bug is added for the genus level group. 

k__kingdom1|p__phylum2|c__class1|o__order1|f__family1|g__genus1|s__species1	5
k__kingdom1|p__phylum2|c__class1|o__order1|f__family1|g__genus1|s__species2	10
add
k__kingdom1|p__phylum2|c__class1|o__order1|f__family1|g__genus1	15

A new kingdom, phylum, class, order, and family is not entered because they would be the same grouping of counts as the new genus level entry.

If we had an additional bug
k__kingdom1|p__phylum2|c__class1|o__order2|f__family12|g__genus23|s__species14	2
k__kingdom1|p__phylum2|c__class1|o__order1|f__family1|g__genus1|s__species1	5
k__kingdom1|p__phylum2|c__class1|o__order1|f__family1|g__genus1|s__species2	10
add
k__kingdom1|p__phylum2|c__class1	17
k__kingdom1|p__phylum2|c__class1|o__order1|f__family1|g__genus1	15

Two new bugs are added because o__order1 and o__order2 can be combined at the c__class1 grouping and s__species1 and s__species2 can be combined at the g__genus1 level. Other groupings at other clade levels are not made because they represent the same groupings of counts already accounted for in the data by bugs and would be redundant. For instance, having a k__kingdom 17 count entry would be the same grouping as the c_class1 bug that was added and so is not created and added.

*How do I reduce the table to a list of bugs?*

> python scripts/scriptManipulateTable.py -i TID -l STSite -b features.txt demo_input/Test.pcl
> python scripts/scriptManipulateTable.py -i TID -l STSite -b 'Bacteria|3417,Bacteria|unclassified|4904' demo_input/Test.pcl

IV. Manipulate samples by metadata.  
J. How do I stratify the table to subtables based on a metadata? (Example. How do I take the HMP table and break it up by body site or time point?)

> python scripts/scriptManipulateTable.py -i TID -l STSite -y STSite demo_input/Test.pcl

K. How do I remove all samples of a certain metadata value? (Example, How do I remove all gut HMP body site samples but leave the rest in the table?)

> python scripts/scriptManipulateTable.py -i TID -l STSite -r STSite,R_Retroauricular_crease, L_Retroauricular_crease demo_input/Test.pcl

V. Manipulate the feature names  
L. How do I add on the 'k__' and 's__' on the names of my bugs?

> python scripts/scriptManipulateTable.py -i TID -l STSite -x demo_input/Test.pcl

VI. Dimensionality Reduction
M. How do I make new composite bugs or metadata using principle components analysis (PCA).
Adding a -p will make principle components (PCs) out of the bug abundance data and then seperately another time with the numeric data in the metadata (This would include boolean indicated as 0 and 1 so make sure to represent them as 'TRUE' and 'FALSE'). The PCA performed is the same as performed by the R function prcomp and was checked by the results of this function. SVD is used in the function and scaling and centering is automatically performed. The feature names are named as such: Metadata/Data + PC# + percent variance which is between 0 and 1 but the decimal was changed to _ so 0.80 would be 0_80. An example would be Data_PC1_0_81234 which means this is a PC made from the abundance data, this was the first PC and the percent variance is 0.81234 or 81% and some change.

> python scripts/scriptManipulateTable.py -i TID -l STSite -p demo_input/Test2.pcl


## scriptPlotFeature.py ##
This script allows you to plot a row of an abundance table, metadata or data. This assumes the first column is the id and the remaining columns are values to be plotted. Three different plots can be generated based on the input arguments and the type of data given to the script. A boxplot is made if two features are given, one numeric and one categorical. A scatterplot is made if two numeric features are given. A histogram is made if one numeric feature is given.  

A. How do I plot a box plot of two data.  
A box plot requires two features, one not categorical and one that is categorical. The script detects this automatically and will plot the correct plot for you as you go.

> python scripts/scriptPlotFeature.py demo_input/Test.pcl STSite 'Bacteria|3417'

B. How do I plot a scatter plot of two data?  
A box plot requires two features, both not categorical. The script detects this automatically and will plot the correct plot for you as you go.

> python scripts/scriptPlotFeature.py demo_input/Test.pcl 'Bacteria|unclassified|4904' 'Bacteria|3417'

C. How do I plot a histogram of a feature?  
Just plot one numeric feature.

> python scripts/scriptPlotFeature.py demo_input/Test.pcl 'Bacteria|3417'

D. How do I change the title or axes?

> python scripts/scriptPlotFeature.py -t Title -x Xaxis -y Yaxis demo_input/Test.pcl 'Bacteria|3417'

E. How do I change the color?  
This script takes hex colors.

> python scripts/scriptPlotFeature.py -c '#333333' demo_input/Test.pcl 'Bacteria|3417'

F. How do i invert the colors for a black background?

> python scripts/scriptPlotFeature.py -r demo_input/Test.pcl 'Bacteria|3417'
