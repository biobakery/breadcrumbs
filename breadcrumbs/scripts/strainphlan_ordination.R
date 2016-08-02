#!/usr/bin/env Rscript

# This script generates an ordination plot from the distance matrix calculated from StrainPhlAn multiple sequence alignment file

# load the optparse library
library(optparse)

# Create command line argument parser
help_description <- "

This script generates an ordination plot from StrainPhlAn output files.

The following positional arguments are required:

clade.distmat_4R.txt (Input file): This is the distmat output file (created by providing the StrainPhlAn MSA file)
metadata.txt (Input file): This is the metadata file
ordination.png (Output file): This is the ordination plot image written"

args <- OptionParser(usage = "%prog clade.distmat_4R.txt metadata.txt ordination.png",
                      add_help_option = TRUE, prog='strainphlan_ordination.R',
                      description=help_description )
args_list <- parse_args(args, positional_arguments=TRUE)

# load libraries after optparse to not load them on help
library(ggplot2)
library(vegan)

# load triangular distance matrix 
e.sir.dist <- read.table( args_list$args[1], sep = "\t", row.names = 1, header = T )
# remove X from colnames
colnames( e.sir.dist ) <- gsub( "X", "", colnames( e.sir.dist ))
# make symmetric, add lower triangle to upper triangle
e.sir.dist[lower.tri(e.sir.dist)] <- t(e.sir.dist)[lower.tri(e.sir.dist)]
# dim(e.sir.dist)  # should be 7 by 7 matrix with demo input files
# ordinate on the distance matrix
e.sir.pcoa <- cmdscale( e.sir.dist, eig = T )
# variance explained 
variance <- head(eigenvals(e.sir.pcoa)/sum(eigenvals(e.sir.pcoa)))
x_variance <- as.integer(variance[1]*100)
y_variance <- as.integer(variance[2]*100)

# get scores for plotting
e.sir.scores <- as.data.frame( e.sir.pcoa$points )

# read in metadata file
e.sir.meta <- read.delim( args_list$args[2], header = T, sep = "\t", row.names = 1 )
# append to e.sir.scores
e.sir.scores.meta <- merge( e.sir.scores, e.sir.meta, by = 'row.names' )
# set first column as rownames and remove it
rownames( e.sir.scores.meta ) <- e.sir.scores.meta[,1]
e.sir.scores.meta[,1] <- NULL
# change colnames
colnames(e.sir.scores.meta) <- c( "PCo1", "PCo2", "SubjectID" )

# plot ordination
png( args_list$args[3], width = 750, height = 600, res = 150 )

ggplot( e.sir.scores.meta, aes(PCo1, PCo2, color=SubjectID) ) + 
  geom_point(size = 4, alpha = 0.75) + theme_classic() + 
  theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank()) + 
  xlab(paste("PCo1 (",x_variance,"% variance explained)")) + ylab(paste("PCo2 (",y_variance,"% variance explained)"))

dev.off()
