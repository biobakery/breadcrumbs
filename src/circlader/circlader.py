#!/usr/bin/env python

#-----------------------------------------------------------------------------
# NAME: circlader.py
# DESCRIPTION:  Circlader (circular cladogram buider) is a python script for
#               creating images of circular cladogram starting from any guide
#               tree in tabular or Newick format
#
# Author: Nicola Segata
# email: nsegata@hsph.harvard.edu
#
# Copyright: (c) 2011
# Licence: <your licence>
#
#-----------------------------------------------------------------------------


import os
import sys,argparse
import circlader_lib as cir

def read_params(args):
    parser = argparse.ArgumentParser(description='Circlader')

    parser.add_argument('tree_file', nargs='?', default=None, type=str,
            help=   "the input tree in Newick format (unless --tf is specified)"
                    "[stdin if not present]")
    parser.add_argument('out_image', nargs='?', default=None, type=str,
            help=   "the output image (the format is guessed from the extension "
                    "[windows visualization if not present]")  
    parser.add_argument('--tree_format', choices=['newick','tabular'], 
                        default='newick', type=str,
            help=       "specifies the input tree format (default \"newick\", "
                        "other choice is \"tabular\")")
    parser.add_argument('--style_file', nargs='?', default=os.getcwd()+"/default_styles/style.txt", type=str,
            help=       "set the style file (default_styles/style.txt if not specified)")
    parser.add_argument('--color_file', nargs='?', default=None, type=str,
            help=       "set the color file (default_styles/colors.txt if not specified)")
    parser.add_argument('--highlight_file', nargs='?', default=None, type=str,
            help=       "set the highlight file (default none)")
    parser.add_argument('--tick_file', nargs='?', default=None, type=str,
            help=       "set the label file for level's names (default none)")
    parser.add_argument('--size_file', nargs='?', default=None, type=str,
            help=       "set the file containing the dimentison of the circles (default none)")
    parser.add_argument('--circle_file', nargs='?', default=None, type=str,
            help=       "set the external circles file (default none) [BETA FEATURE]")
    parser.add_argument('--format', choices=['png','pdf','ps','eps','svg'], default=None, type=str,
            help=       "set the format of the output image (default none "
                        "meaning that the format is guessed from the output "
                        "file extension)")
    parser.add_argument('--dpi', default=300, type=int )

    return vars(parser.parse_args())

params = read_params(sys.argv)

cladogram = cir.Tree()
cladogram.read_colors(params['color_file'])
cladogram.read_style(params['style_file'])
cladogram.read_sizes(params['size_file'])
cladogram.read_circles(params['circle_file'])
cladogram.read_highlights(params['highlight_file'])
cladogram.read_tick_labels(params['tick_file'])
if params['tree_format'] == 'newick':
    cladogram.load_newick(params['tree_file'])
else:
    cladogram.load_lefse(params['tree_file'])
cladogram.pos_rad_leaves()
cladogram.set_pos()
cladogram.draw(params['out_image'],outformat=params['format'],dpi=params['dpi'])
