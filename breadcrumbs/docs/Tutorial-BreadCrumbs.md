# BreadCrumbs Tutorial #

This is a brief tutorial to get you acquainted with the scripts provided in breadcrumbs. This tutorial is oragnized by script and task. Examples are given using files in the demo_input folder which is included in the BreadCrumbs package. Each of these commands should work from the command line in the breadcrumbs directory.

Please note all of the following calls expect you to be in the breadcrumbs directory (the relative path to the demo_input files will need to be modified otherwise).

Enjoy and happy science!

## Contents: ##
1. scriptPCoA  
2. scriptManipulateTable.py  
I. Manipulating the measurements  
II. Filtering  
III. Filtering with knowledge of feature hierarchical relationship
IV. Manipulate samples by metadata
V. Manipulate the feature names
3. scriptPlotFeature.py
4. scriptBiplotTSV.R
5. scriptConvertBetweenBIOMAndPCL.py

## scriptPCoA.py ##
This script allows one to plot a PCoA of an abundance table. In the plot each sample is one marker. The marker shape and color is determined by a metadata (of your choice). The distances between each sample is determined by a specific beta-diversity distance metric. By default Bray-curtis distance is used. This can be changed as needed. You will notice for every call you must give it the sample id (-i) and the last metadata which should be the row before your first data (-l). This helps the scripts understand what is a data measurement and what is a metadata.

A. How do I make a PCoA of an abundance table, painting (coloring) it by a specific metadata?

> scriptPcoa.py -i TID -l STSite -p STSite demo_input/Test.pcl

B. How do I make a series of PCoAs of an abundance table, one PCoA for every metadata?

If nothing is specified with -p then all metadata are painted. Note there are a max of 9 shapes to use, a metadata will be skipped if it has more than 9 levels (specific values which can be used many times). Don't worry, the script will let you know if this happens and will just skip to the next metadata.

> scriptPcoa.py -i TID -l STSite demo_input/Test.pcl

C. How do I use a different beta-diversity distance metric instead of Bray-curtis distance?
The following metrics can be choosen: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, euclidean, hamming, sqeuclidean, unifrac_unweighted, unifrac_weighted

> scriptPcoa.py -i TID -l STSite -m sqeuclidean demo_input/Test.pcl

D. How do I get the coordinates of the points in the PCoA plot? Use -C and give a file path to which to write.

> scriptPcoa.py -i TID -l STSite -C coordinates.txt demo_input/Test.pcl

E. How do I get the distance matrix represented by the PCoA plot? Use -D and give a file path to which to write.

> scriptPcoa.py -i TID -l STSite -D distances.txt demo_input/Test.pcl

F. How do I make a PCoA using unifrac type metrics.

> scriptPcoa.py -m unifrac_weighted -t demo_input/GreenGenesCore-May09.ref.tre -e demo_input/fastunifrac_Ley_et_al_NRM_2_sample_id_map.txt -c demo_input/fastunifrac_Ley_et_al_NRM_2_sample_id_map-colors.txt
> scriptPcoa.py -m unifrac_unweighted -t demo_input/GreenGenesCore-May09.ref.tre -e demo_input/fastunifrac_Ley_et_al_NRM_2_sample_id_map.txt -c demo_input/fastunifrac_Ley_et_al_NRM_2_sample_id_map-colors.txt

G. How do I change the marker transparency?
Use -a , --alpha . Change alpha to .5 (half) remember 0-1 are the range of values, 0 = invisible, 1 = no transparency

> scriptPcoa.py -a .5 -i TID -l STSite -p STSite -o TestChangeAlpha.pdf demo_input/Test.pcl

H. How do I turn off marker outlining (black outlining line)
Use -u or --outlining
> scriptPcoa.py -u -i TID -l STSite -p STSite -o TestNoOutline.pdf demo_input/Test.pcl

I. How do i make the markers bigger in size? Use the -z or --size argument with an integer value (default 20). Bigger numbers are bigger shapes.

> scriptPcoa.py -z 40 -i TID -l STSite -p STSite -o TestSizeBigMarkers.pdf demo_input/Test.pcl

J. How do I make the markers smaller in size? Use the -z or --size argument with an integer value (default 20). Smaller numbers are smaller shapes, no negative numbers please.

> scriptPcoa.py -z 5 -i TID -l STSite -p STSite -o TestSizeSmallMarkers.pdf demo_input/Test.pcl

K. How do I move the lengend?
Use the -g or --legend (valid locations are upper_right, upper_left, lower_left, lower_right, right, center_left, center_right, lower_center, upper_center, center)

> scriptPcoa.py -g lower_right -i TID -l STSite -p STSite -o TestLegend.pdf demo_input/Test.pcl

L. How do I turn off the shapes?
Use -S, --noShape to make all the markers circles.

> scriptPcoa.py -S -i TID -l STSite -p STSite -o TestNoShape.pdf demo_input/Test.pcl

M. How do I force the coloring of the markers?
Use the -r, --customColor and give colors based on the legend labels. Pairs are label:color and are comma delimited, no white spaces please!

> scriptPcoa.py -r "R_Retroauricular_crease:red,R_Antecubital_fossa:orange,L_Antecubital_fossa:yellow,L_Retroauricular_crease:green,Subgingival_plaque:blue" -i TID -l STSite -p STSite -o TestCustomColor.pdf demo_input/Test.pcl


There already exists a collection of functionality surrounding unifrac distances in Qiime and related software. We support these metrics here for completeness, if your need is not met here, please look into Qiime and related software for a solutions with a more rich collection of functionality.

## scriptManipulateTable.py ##
Abundance tables can be difficult to manipulate. This script captures frequent tasks that may be important to manipulating an abundance table including normalization, summing, filtering, stratifying the tables into subsets (for instance breaking up a large HMP table into tables, one for each body site), and other functionality. You will notice for every call you must give it the sample id (-i) and the last metadata which should be the row before your first data (-l). This helps the scripts understand what is a data measurement and what is a metadata.

_Remember you can do multiple tasks or use multiple arguments at the same time._
_Here is an example of summing, normalizing, adding on clade prefixes, and stratifies the tables based on the STSite metadata_

> scriptManipulateTable.py -i TID -l STSite -s -n -x -y STSite demo_input/Test.pcl

Please look at the detailed description of normalization and summation for a clear understanding of how the data is being manipulated.

*Manipulating the measurements*  
A. How do I sum a table based on clade names?

> scriptManipulateTable.py -i TID -l STSite -s demo_input/Test.pcl

B. How do I normalize a table?

> scriptManipulateTable.py -i TID -l STSite -n demo_input/Test.pcl

*Filtering*  
C. How do I filter a normalized table by percentage?

This filters out bugs that are not in the top 95 percentage of at least 5 percent of the samples (a good default).

> scriptManipulateTable.py -i TID -l STSite -P 95,5 demo_input/Test.pcl

D. How do I filter a normalized table by a minimum abundance?  

This filters out bugs that do not have at least a certain number of bugs in a certain number of samples. Here we show
the call to filter out all bugs which do not have at least 3 samples with at least 0.0001 abundance (a good initial default).

> scriptManipulateTable.py -i TID -l STSite -A 0.0001,3 demo_input/Test.pcl

E. How do I filter a count table by count occurrence?  

This removes samples that do not have at least 5 counts in at least 3 samples (an initial default to use could be 2,2).

> scriptManipulateTable.py -i TID -l STSite -O 5,3 demo_input/Test.pcl

F. How do I filter a table by standard deviation?  

> scriptManipulateTable.py -i TID -l STSite -D 1 demo_input/Test.pcl

*Filtering with knowledge of feature hierarchical relationship.*  
F. How do I make the table have only terminal nodes?

> scriptManipulateTable.py -i TID -l STSite -t demo_input/Test.pcl

G. How do I remove all the OTUs from a table?

> scriptManipulateTable.py -i TID -l STSite -u demo_input/Test.pcl

H. How do I reduce all bugs more specific than a certain clade? Aka, how do I reset a table to be only a clade (genus) or higher?  

This reduces the bugs to bugs with 3 levels of hierarchy or less (class on a standard biological taxonomy).

> scriptManipulateTable.py -i TID -l STSite -c 3 demo_input/Test.pcl

You may want to hierarchically sum all of you bugs before reducing the table to a certain level, just in case you are missing some.

> scriptManipulateTable.py -i TID -l STSite -s -c 3 demo_input/Test.pcl

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

> scriptManipulateTable.py -i TID -l STSite -b features.txt demo_input/Test.pcl
> scriptManipulateTable.py -i TID -l STSite -b 'Bacteria|3417,Bacteria|unclassified|4904' demo_input/Test.pcl

IV. Manipulate samples by metadata.  
J. How do I stratify the table to subtables based on a metadata? (Example. How do I take the HMP table and break it up by body site or time point?)

> scriptManipulateTable.py -i TID -l STSite -y STSite demo_input/Test.pcl

K. How do I remove all samples of a certain metadata value? (Example, How do I remove all gut HMP body site samples but leave the rest in the table?)

> scriptManipulateTable.py -i TID -l STSite -r STSite,R_Retroauricular_crease, L_Retroauricular_crease demo_input/Test.pcl

V. Manipulate the feature names  
L. How do I add on the 'k__' and 's__' on the names of my bugs?

> scriptManipulateTable.py -i TID -l STSite -x demo_input/Test.pcl

## scriptPlotFeature.py ##
This script allows you to plot a row of an abundance table, metadata or data. This assumes the first column is the id and the remaining columns are values to be plotted. Three different plots can be generated based on the input arguments and the type of data given to the script. A boxplot is made if two features are given, one numeric and one categorical. A scatterplot is made if two numeric features are given. A histogram is made if one numeric feature is given.  

A. How do I plot a box plot of two data.  
A box plot requires two features, one not categorical and one that is categorical. The script detects this automatically and will plot the correct plot for you as you go.

> scriptPlotFeature.py demo_input/Test.pcl STSite 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|72'

B. How do I plot a scatter plot of two data?  
A box plot requires two features, both not categorical. The script detects this automatically and will plot the correct plot for you as you go.

> scriptPlotFeature.py demo_input/Test.pcl 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|4904' 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|72'

C. How do I plot a histogram of a feature?  
Just plot one numeric feature.

> scriptPlotFeature.py demo_input/Test.pcl 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|72'

D. How do I change the title or axes?

> scriptPlotFeature.py -t Title -x Xaxis -y Yaxis demo_input/Test.pcl 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|72'

E. How do I change the color?  Use -c and a hex color.

> scriptPlotFeature.py -c '#333333' demo_input/Test.pcl 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|72'

F. How do I invert the colors for a black background? Use -r .

> scriptPlotFeature.py -r demo_input/Test.pcl 'Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|unclassified|72'

## scriptBiplotTSV.R ##
This script allows one to plot a tsv file as a biplot. A tsv file is a transposed PCL file (demo files are found in demo_input). The positioning of sample markers and bug text are generated by nonmetric multidimensional scaling (R Vegan package). The metadata are represented by arrows and then a text at the head of the arrow. The coordinates of the arrows are determined by the center/average of the coordinates of the samples with that metadata showing a central tendency of where that metadata is located. More specifically, discontinuous metadata are broken down to levels (values), then each level is made into it's own binary metadata (0 for not having that value and 1 for having that value). Then for each new metadata, samples with the value of 1 are selected and have their coordinates in the ordination are averaged. This average coordinate set is then used as the coordinates for that metadata level. For continuous data, using the ordination coordinates for all the sample points, the value of the continuous metadata is placed in a landscape using the sample coordinates as x and y and the z as the metadata value. This is then smoothed with a LOESS and then the maximum fitted value's coordinates are used as the central tendency of the metadata.

This is a customizable plot there the metadata plotted, the bugs plotted, the arrow color, the bug and metadata text color, the sample colors, the sample shapes, and the title can be changed. Below are examples of how to use the commandline for this script; although options are shown seperately here, many options can be used together.

A. This is the minimal call. The call to the script must include the following positional arguments lastmetadata and inputPCLFile after any optional arguments (given with flags preceeded by - or --).

> scriptBiplotTSV.R STSite demo_input/Test-Biplot.tsv

This consists of the script call, the last metadata value, and the input file.

B. How do I specify an output file name? Use -o and then a name ending with the extension .pdf .

> scriptBiplotTSV.R -o Test2Biplot.pdf STSite demo_input/Test-Biplot.tsv

C. How do I specify bug names to plot? Use -b and the names of the bugs to plot (as written in your pcl file) seperated by commas.

> scriptBiplotTSV.R -b 'Bacteria|3417' STSite demo_input/Test-Biplot.tsv
> scriptBiplotTSV.R -b 'Bacteria|3417,Bacteria|unclassified|4904' STSite demo_input/Test-Biplot.tsv

D. How do I specify metadata to plot? Using the -m option, if plotting a continuous metadata, give the id of the metadata (first column entry of metadata), for any other metadata concatonate the metadata ID and the value of interest with "_" . Here are a working examples.

> scriptBiplotTSV.R -m 'STSite_L_Antecubital_fossa,STSite_R_Antecubital_fossa' STSite demo_input/Test-Biplot.tsv
> scriptBiplotTSV.R -m 'Continuous' STSite demo_input/Test-Biplot.tsv

E. How do I specify a title? Use -i and the title text. Use ' to surround the text. This helps the command line understand that the text is together and not a series of commands. These should be used when you are giving a flag a value that has spaces or anything but alphanumeric characters.

> scriptBiplotTSV.R -i 'Test Title' STSite demo_input/Test-Biplot.tsv

F. How do I specify a metadata to shape markers by? Use -y and the id for your metadata in your pcl file (the entry for your metadata in the first column).

> scriptBiplotTSV.R -y STSite STSite demo_input/Test-Biplot.tsv

G. How do I specify specific shapes to use? This requires the combination of -y and -s . Use -y to specify the metadata to use. Use -s to specify what shapes should be used for what metadata values. These are given as metadatavalue:shape,metadataValue:shape

> scriptBiplotTSV.R -y STSite -s 'L_Antecubital_fossa:15,R_Antecubital_fossa:23' STSite demo_input/Test-Biplot.tsv

H. How do I specify a metadata to color markers by? Use -c and the id for your metadata in your pcl file (the entry for your metadata in the first column).

> scriptBiplotTSV.R -c STSite STSite demo_input/Test-Biplot.tsv

I. How do I specify a default marker shape to use instead of using a metadata? Use -d and a number recognized by R's pch parameter (number between 1-25). For more information http://www.statmethods.net/advgraphs/parameters.html

> scriptBiplotTSV.R -d 1 STSite demo_input/Test-Biplot.tsv

J. How do I specify a color range to use when coloring? Use -r with two (R supported) colors seperated by a comma. R supported colors can be found in many sources including this one http://www.stats4stem.org/r-colors.html

> scriptBiplotTSV.R -r 'red,cyan' STSite demo_input/Test-Biplot.tsv

K. How do I specify a color to use when drawing arrows? Use -a and a (R supported) color. R supported colors can be found in many sources including this one http://www.stats4stem.org/r-colors.html

> scriptBiplotTSV.R -a  red STSite demo_input/Test-Biplot.tsv

L. How do I specify a color to use for arrow text? Use -w and a (R supported) color. R supported colors can be found in many sources including this one http://www.stats4stem.org/r-colors.html

> scriptBiplotTSV.R -w orange STSite demo_input/Test-Biplot.tsv

M. How do I specify a color to use for bug text? Use -t and a (R supported) color. Make sure to use -b to plot bugs.  R supported colors can be found in many sources including this one http://www.stats4stem.org/r-colors.html

> scriptBiplotTSV.R -t pink -b 'Bacteria|3417' STSite demo_input/Test-Biplot.tsv

N. How do I rotate the projection (plot) in reference to a specific metadata? Use the -e option and give the plot a metadata and a weight for that metadata, the larger the weight, the more the rotation takes into account the metadata. You may have to experiement with different weights and see how the rotation is affected. The weights can be from 0 (no rotation by the metadata) to a very large number. The metadata name should be the metadata id if the value is continuous or the metadata id and the value (level) of interest seperated by a _ if the metadata is not continuous. Below are two examples, one using a continuous metadata and one using a discontinuous metadata.

> scriptBiplotTSV.R -e 'Continuous,2' STSite demo_input/Test-Biplot.tsv
> scriptBiplotTSV.R -e 'STSite_L_Antecubital_fossa,.5' STSite demo_input/Test-Biplot.tsv

O. How do I color NAs a specific color, no matter other coloring in the plot? Use -n and a color supported by R.  R supported colors can be found in many sources including this one http://www.stats4stem.org/r-colors.html This requires you to be coloring the plot by a metadata (option -c).

> scriptBiplotTSV.R -n grey -c STSite STSite demo_input/Test-BiplotNA.tsv

P. How do I scale arrows in the plot. Use -z and a number to weight how much the metadata influences the rotation (number between 0 and very large).

> scriptBiplotTSV.R -z 2 STSite demo_input/Test-Biplot.tsv

Q. How do I plot metadata labels without the arrows?

> scriptBiplotTSV.R -A STSite demo_input/Test-Biplot.tsv

R. How do I plot the biplot without metadata?

> scriptBiplotTSV.R -m "" STSite demo_input/Test-Biplot.tsv

## scriptConvertBetweenBIOMAndPCL.py ##
The script allows one to convert between PCL and BIOM file formats. ID, last feature (row) metadata, and last sample metadata are optional information in the script call (when converting from PCL to BIOM). These are used to dictate placement of certain key sample metadata in the PCL file. Typically, it is helpful to set these arguments. This aids in the consistent and reliable manipulation of these files. If the are not given, a guess will be made to the ID and it will be assumed no metadata exist.

A quick definition:
*ID or sample id* -  typically your first row in the PCL file (the Ids of all your samples) in the example below "ID"
*Feature (row) metadata* - columns in your PCL file which describe your features. These come after your feature IDs but before your measurements.
*Sample metadata* - rows in your PCL file which come before your measurements and describe your samples

For a description of a PCL and it's parts please look in the docs folder for PCL-Description.txt

A. The minimal call to convert from BIOM file to a PCL file or visa versa. This call indicates the sample metadata entry which is the sample id and which is the last listed metadata in a pcl file (before the data measurements). When converting a PCL file, if there are no metadata and only a metadata id, -l  and -i is not required. If there are multiple metadata in a pcl file the -l (last metadata) field is required. Neither of these fields are required for biom file conversion to pcl.

> scriptConvertBetweenBIOMAndPCL.py demo_input/Test_no_metadata.pcl example1.biom
> scriptConvertBetweenBIOMAndPCL.py demo_input/Test.biom example2.pcl
> scriptConvertBetweenBIOMAndPCL.py -l STSite demo_input/Test.pcl example3.biom

B. Specifying ID and lastmetadata

> scriptConvertBetweenBIOMAndPCL.py -i TID -l STSite demo_input/Test.pcl example4.biom
> scriptConvertBetweenBIOMAndPCL.py -i TID -l STSite demo_input/Test.biom example5.pcl

C. The case where there are no sample metadata, just sample IDs. Indicate the ID and if no last metadata is indicated (-l) it is assumed no sample metadata exist.

> scriptConvertBetweenBIOMAndPCL.py -i ID demo_input/Test_no_metadata.pcl example6.biom
> scriptConvertBetweenBIOMAndPCL.py -i ID demo_input/Test_no_metadata.biom example7.pcl

D. The case when converting a PCL file with Feature (row) metadata (for example taxonomy_5). Include the last column with feature metadata.

> scriptConvertBetweenBIOMAndPCL.py -i ID -r taxonomy_5 -l STSite demo_input/testFeatureMetadata.pcl testFeatureMetadata.biom

E. Although the output file name can be automatically generated, the output file name can be given if needed.

> scriptConvertBetweenBIOMAndPCL.py -i TID -l STSite demo_input/Test.biom CustomFileName.pcl
> scriptConvertBetweenBIOMAndPCL.py -i TID -l STSite demo_input/Test.pcl CustomFileName.biom

F. Indicate the use of a pcl file using a delimiter that is not tab or indicate the creation of a pcl file using a delimier that is not tab.

> scriptConvertBetweenBIOMAndPCL.py -i TID -l STSite -f , demo_input/Test-comma.pcl
> scriptConvertBetweenBIOMAndPCL.py -i TID -l STSite -f , demo_input/Test-comma.biom

G. Piping the input and output data. If piping is used (-p) one must pipe in and out in stdin and stdout. As well, indicate what format the file should be written in with -m . Here is an example

> cat demo_input/Test_no_metadata.pcl | scriptConvertBetweenBIOMAndPCL.py -p -m biom
> cat demo_input/Test_no_metadata.biom | scriptConvertBetweenBIOMAndPCL.py -p -m pcl
