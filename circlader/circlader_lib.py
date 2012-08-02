#-----------------------------------------------------------------------------
# NAME: circlader_lib.py
# DESCRIPTION:  Circlader (CIRcular CLADogram buidER) is a python script for
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

import sys,os,math,matplotlib
#matplotlib.use('TkAgg')                                         
matplotlib.use('Agg')                                         
#matplotlib.use('PDF')
from matplotlib import collections
from Bio import Phylo
import numpy as np
from pylab import *
import operator
import matplotlib.patches as mpatches

class Tree:
    max_rad_dist = math.pi*0.15

# Class specifying clade with strutural characteristics and 
# associated information
    class Clade:
        def __init__(   self,taxa_id=0, name='',
                        br_len=-1.0, root_br_dist=-1.0,
                        highlighted=False,
                        ext_seg = False,
                        ext_seg_vec = None):
            self.id = taxa_id
            self.name = name
            self.label = name
            self.label_cat = ''
            self.col = 'w'
            self.br_len = br_len
            self.tot_br_len = 0.0 
            self.root_br_dist = root_br_dist
            self.is_leaf = True
            self.is_highlighted = highlighted
            self.__children = {}
            self.pos = Tree.VisPos()
            self.nleaves = 0
            self.size = 0
            self.ext_seg = ext_seg
            self.ext_seg_vec = ext_seg_vec

        def add_child(self,cl):
            self.__children[cl.id] = cl
            self.tot_br_len += (cl.br_len + 
                    (0 if cl.is_leaf else cl.tot_br_len))
            self.is_leaf = False
            
        def get_children(self):
            return self.__children.values()

#   Class decribing graphical information associated with clades
    class VisPos:
        def __init__(   self, r=0.0, rad=0.0, rmin=0.0, 
                        rmax=0.0, lab_lev = 0.0):
            self.rad = rad
            self.r = r
            self.rad_min = rmin
            self.rad_max = rmax
            self.lab_lev = lab_lev
    def __init__(self):
        self.__all_taxa = set([])
        self.__noname_gen = self.__unique_none_id()
        self.__max_br_depth = 0.0
        self.__min_br_depth = float(sys.maxint) 
        self.__leaves = []
        self.min_high = float(sys.maxint)
        self.max_high = -1.0
        self.added_label_cats = []
        self.vis_cats = []
        self.wing_ext_max = 1.0
        self.cseg = None 

    def __unique_none_id(self):
        for i in xrange(2**31-1):
            yield "_"+str(i)


    def add_clade(self,cl,fn,br_depth):
        cl.size = self.sizes[cl.id] if hasattr(self, 'sizes') and cl.id in self.sizes else self.opt['default_taxa_size']
        if cl.is_highlighted:
            cl.pos.lab_lev = br_depth+1.0
            if self.min_high > cl.pos.lab_lev:
                self.min_high = cl.pos.lab_lev
            if self.max_high < min( cl.pos.lab_lev,
                    float(self.opt['highlighting_bar_stop_level'])*1.001):
                self.max_high = cl.pos.lab_lev
            cl.label = self.labels[fn]
            cl.col = self.label_color[fn]
            cl.label_cat = self.label_cat[fn]
            cl.size *= self.opt['highlight_taxa_size_magnifier']


    def load_lefse(self,inp_f):
        with open(inp_f, 'r') as inp:
            rows = ["root."+l.rstrip().split("\t")[0] for l in inp.readlines()]
        rows.append("root")

        self.opt['ignore_branch_len'] = 1

        def rec_add(clade_name,rows,br_depth=0.0, first = False):
            self.__all_taxa.add(clade_name)
            fn = clade_name if not clade_name.startswith("root.") \
                            else clade_name.split("root.")[1]
            highlighted = (fn in self.labels)
            ext_seg = (self.cseg and fn in self.cseg)
            ext_seg_v = self.cseg[fn] if ext_seg else None
            cl = Tree.Clade(    clade_name, clade_name.split(".")[-1],1.0,
                                br_depth+1.0, highlighted = highlighted, ext_seg = ext_seg, ext_seg_vec = ext_seg_v )
            self.add_clade(cl,fn,br_depth)
            rows.remove(clade_name)
            np = clade_name.count(".")
            children = [r for r in rows if r.count(".") == np+1 and \
                        r.count(".") > 0  and \
                        ".".join(r.split(".")[:-1]) == clade_name]
            if not children:
                self.__leaves.append(cl)
                if cl.root_br_dist > self.__max_br_depth:
                    self.__max_br_depth = cl.root_br_dist
                if cl.root_br_dist > self.__max_br_depth:
                    self.__max_br_depth = cl.root_br_dist
                cl.nleaves = 1
            nleav = len(self.__leaves)
            for c in children:
                cl.add_child(rec_add(c,rows,br_depth+1.0 if not first else 0.0))
            sep = 'sep_clades' in self.opt and self.opt['sep_clades']
            if sep and children and all([c.is_leaf for c in cl.get_children()]):
                self.__leaves.insert(nleav,None)
                self.__leaves.append(None)
            cl.nleaves = sum([c.nleaves for c in cl.get_children()])

            return cl

        self.root = rec_add(rows[-1],rows, first = True)    

    def load_newick(self,inp_f):
        if os.path.splitext(inp_f)[-1] == '.nwk':
            bio_tree = Phylo.read(inp_f, "newick")
        elif os.path.splitext(inp_f)[-1] == '.xml':
            bio_tree = Phylo.read(inp_f, "phyloxml")
        else:
            sys.stderr.write( "Unrecognized tree extensions: "+os.path.splitext(inp_f)[-1]+"\n" )
            sys.exit(0)

        def rec_add(clade,br_depth=-0.0):
            nam = clade.name if clade.name else "noname"+self.__noname_gen.next()
            if nam in self.__all_taxa:
                oldn = nam
                nam = nam+self.__noname_gen.next()
                print "Warning: "+oldn+" non unique, renamed to "+nam
            self.__all_taxa.add(nam)
            if self.opt['ignore_branch_len']: clade.branch_length = 1.0
            fn = nam if not nam.startswith("root.") else nam.split("root.")[1]
            highlighted = (fn in self.labels)
            ext_seg = (self.cseg and fn in self.cseg)
            ext_seg_v = self.cseg[fn] if ext_seg else None
            cl = Tree.Clade(    nam,
                                nam,
                                clade.branch_length,
                                br_depth+clade.branch_length,
                                highlighted=highlighted,
                                ext_seg = ext_seg, ext_seg_vec = ext_seg_v)
            self.add_clade(cl,fn,br_depth)
            
            if not clade.clades: 
                cl.nleaves = 1
                return cl 

            if br_depth+cl.br_len < self.opt['max_branch_depth']:
                nleav = len(clade.clades) 
                for c in clade.clades:
                    cl.add_child(rec_add(c,br_depth+cl.br_len))
            cl.nleaves = sum([c.nleaves for c in cl.get_children()])
            return cl

        def rec_leaves( cl ):
            children = cl.get_children()
            if 'size_sorted' in self.opt and self.opt['size_sorted']:
                children = sorted(children, key=lambda x: x.nleaves,reverse= True)
            
            for c in children:
                rec_leaves( c ) 

            if not children:
                self.__leaves.append(cl)
                if cl.root_br_dist > self.__max_br_depth:
                    self.__max_br_depth = cl.root_br_dist
                if cl.root_br_dist < self.__min_br_depth:
                    self.__min_br_depth = cl.root_br_dist
                cl.nleaves = 1
                return
            nleav = len(self.__leaves)
            sep = 'sep_clades' in self.opt and self.opt['sep_clades'] 
            
            if sep and all([c.is_leaf for c in children]):
                self.__leaves.insert(nleav,None)
                self.__leaves.append(None)

        self.root = rec_add(bio_tree.root,-1.0)

        rec_leaves( self.root )
    
    def pos_rad_leaves(self):
        nl = len(self.__leaves)
        cir_width = float(self.opt['total_rotation'])/360.0 
        self.cir_offset = np.pi*0.5+(1.0-cir_width)*np.pi
        dist = math.pi*cir_width*2.0/nl
        if dist > self.__class__.max_rad_dist:
            dist = self.__class__.max_rad_dist
        prev = self.cir_offset-dist
        ln = False
        for i,c in enumerate(self.__leaves):
            if not c:
                ln = True
                continue
            c.pos.rad = dist*(i if ln else i)+self.cir_offset
            c.pos.rad_min = (c.pos.rad + prev)*0.5 
            next_none = self.__leaves[(i+1)%len(self.__leaves)]
            c.pos.rad_max = (c.pos.rad + dist*(i+1 if next_none else i+3) +self.cir_offset)*0.5  
            prev = c.pos.rad
            ln = False
    
    def set_pos(self):
        cl = self.root
        def rec_set_pos(clade):
            children = clade.get_children()
            #children = sorted(children, lambda x,y: cmp(x.nleaves,y.nleaves),reverse= True)
            clade.pos.r = clade.root_br_dist / self.__max_br_depth
            if children:
                rad_mm = [rec_set_pos(c) for c in children]
                rads = [lcl.pos.rad for lcl in children]
                clade.pos.rad = (min(rads)+max(rads))*0.5
                clade.pos.rad_min,clade.pos.rad_max = min([r[0] for r in rad_mm]), max([r[1] for r in rad_mm]) 
                return clade.pos.rad_min,clade.pos.rad_max
            return clade.pos.rad_min,clade.pos.rad_max
        rec_set_pos(cl)


    def draw(self,out_img_file,outformat=None,dpi=300):
        fig = plt.figure(figsize=(7,7))
        ax = fig.add_subplot(   111, 
                                polar=True, 
                                frame_on=False )

        plt.subplots_adjust(    right=1.0-self.opt['right_space'],
                                left=self.opt['left_space'],
                                top=1.0-self.opt['top_space'],
                                bottom=self.opt['bottom_space'] )
        xticks([])
        yticks([])

        circle_taxa_draw = []

        ign_bl = self.opt['ignore_branch_len']
        ext_unit = self.opt['annotation_wing_sep'] 
        high_start = float(self.opt['highlighting_bar_start_level'])/self.__max_br_depth
        high_stop = float(self.opt['highlighting_bar_stop_level']+1)/self.__max_br_depth*1.001

        taxa_circle, taxa_circle_h = {}, {}
        taxa_circle_attr = [ 'rad', 'r', 'shape', 'col', 'size', 'linew' ]
        for a in taxa_circle_attr:
            taxa_circle[a], taxa_circle_h[a] = [], []

        branch_line_1, branch_line_2, branch_line_3 = {}, {}, {}
        branch_line_attr = [ 'frto', 'col' ]
        for a in branch_line_attr:
            branch_line_1[a], branch_line_2[a], branch_line_3[a]  = [], [], []

        def draw_branch(clade,children,fcol):
            rads = [lcl.pos.rad for lcl in children]
            min_rads,max_rads = min(rads),max(rads) 
            sbl = self.opt['sub_branches_opening_point'] if ign_bl else 1.0
            sb,rsb = sbl,1-sbl 
            redf = 1.0-self.opt['sub_branches_angle_reduction'] 
            red,nred = (redf,1-redf) if abs(max_rads - min_rads) < math.pi else (1.0,0.0)

            mid = ( max_rads + min_rads ) * 0.5 
        
            rads_l = list(arange( min_rads*red+mid*nred, max_rads*red+mid*nred,0.05)) \
                        + [max_rads*red+mid*nred]

            
            if self.opt['highligh_branch_color']:
                col = fcol 
                if clade.is_highlighted:
                    col = clade.col
            else: col = self.opt['branch_color']

            if clade != self.root:
                branch_line_1['frto'].append(   
                        np.array( 
                            [np.array(
                                [c, sb*clade.pos.r+rsb*children[0].pos.r] ) 
                                for c in rads_l]   )   )
                branch_line_1['col'].append(col)
            
            branch_line_2['frto'].append(   
                    np.array( 
                        [np.array( [clade.pos.rad*red+mid*nred, sb*clade.pos.r+rsb*children[0].pos.r]),
                        np.array( [clade.pos.rad, clade.pos.r])] ) )  

            branch_line_2['col'].append(col)
            
            """
            if clade == self.root:
                for c in children:
                    branch_line_3['col'].append(col)
                    print [c.pos.rad, c.pos.r]
                    branch_line_3['frto'].append( np.array( [   np.array( [c.pos.rad*red+mid*nred, sb*clade.pos.r+rsb*c.pos.r] ),
                                                  np.array( [c.pos.rad, c.pos.r] ) ] ) ) 
#branch_line_3['frto'].append( np.array( [   np.array( [c.pos.rad, c.pos.r] ), # it should be [0.0,0.0] but with that the branches disappear
#                                                                np.array( [c.pos.rad, 0.0] ) ] ) )
            """  
            if clade == self.root:
                for c in children:
                    branch_line_3['col'].append(col)
                    branch_line_3['frto'].append( np.array( [   np.array( [c.pos.rad*red+mid*nred, 0.0] ),
                                                                np.array( [c.pos.rad, c.pos.r] ) ] ) )
            else:
                for c in children:
                    branch_line_3['col'].append(col)
                    branch_line_3['frto'].append( np.array( [   np.array( [c.pos.rad*red+mid*nred, sb*clade.pos.r+rsb*c.pos.r] ),
                                                            np.array( [c.pos.rad, c.pos.r] ) ] ) )
            return col
            

        def draw_taxa_circle(clade):
            if clade.is_leaf and not self.opt['draw_leaf_taxa']:
                return
            if not clade.is_leaf and not self.opt['draw_internal_taxa']:
                return
            tc = taxa_circle_h if clade.is_highlighted else taxa_circle
            tc['rad'].append( clade.pos.rad )
            tc['r'].append( clade.pos.r )
            tc['shape'].append( 'o' )
            tc['col'].append( clade.col )
            tc['size'].append( clade.size  )
            tc['linew'].append( self.opt['taxa_circle_edge_width'] )

        def draw_extended_leaves(clade):
            ax.plot(    [clade.pos.rad,clade.pos.rad],
                        [clade.pos.r,1.0],"-",color=[0.7,0.7,0.7])
 
        def draw_wing(clade):
            ext = self.opt['fixed_wing_stop'] if 'fixed_wing_stop' in self.opt and self.opt['fixed_wing_stop'] >= 0 else self.max_high-clade.pos.lab_lev+1
            ext_max = ext*ext_unit
            
            if not clade.label_cat in [v['label'] for v in self.vis_cats] and clade.label_cat:
                self.vis_cats.append(   {   'color':clade.col,
                                            'alpha':self.opt['bar_alpha'],
                                            'label':clade.label_cat}    )
            
            if not clade.label and 'no_wind_if_no_label' in self.opt and self.opt['no_wind_if_no_label']:
                return

            if high_start < clade.pos.r <= high_stop:
                wlab = clade.label.split("_r_")[-1] if not clade.label.count(":") else clade.label.split(":")[0].split("_r_")[-1]

                ax.bar(     clade.pos.rad_min, 
                            1-clade.pos.r+ext_max+ext_unit*0.25,
                            width = abs(clade.pos.rad_min-clade.pos.rad_max),
                            bottom = clade.pos.r,
                            alpha = self.opt['bar_alpha'],
                            color=clade.col, 
                            edgecolor=clade.col )
                if 1.0+ext_max+ext_unit*0.25 > self.wing_ext_max:
                    self.wing_ext_max = 1.0+ext_max+ext_unit*0.25
                

                des = float(180.0*(clade.pos.rad_min+clade.pos.rad_max)/np.pi)*0.5 # -(0 if clade.label.count("_r_") else 90)
                if not clade.label.count("_r_"):
                    des += (90 if 180 < des < 360 else -90)

                lro = 0.5 if 'radial_label_wing_offset' not in self.opt else self.opt['radial_label_wing_offset']
                
                ax.text(    ((clade.pos.rad_min+clade.pos.rad_max)*0.5), 
                            1+ext_max-ext_unit*lro+ext_unit*0.25,
                            wlab,
                            rotation=des,
                            ha="center",
                            va="center",
                            fontsize=self.opt['label_font_size'], 
                            fontstretch=0   )
                if clade.label.count(":"):
                    ax.bar( 0.0, 
                            0.0, 
                            width = 0.0, 
                            bottom = 0.0, 
                            alpha = 1.0, 
                            color=clade.col, 
                            label=clade.label.split("_r_")[-1]   )


        def draw_sectors(clade):
            for i,v in enumerate(clade.ext_seg_vec):
                if v[1] <= 0.0: continue
                height = self.opt['seg_radial_depth'] * ( v[3] if v[3] else 1.0)
                width = abs(clade.pos.rad_min-clade.pos.rad_max)*self.opt['seg_width']
                startx = clade.pos.rad_min+abs(clade.pos.rad_min-clade.pos.rad_max)*((1.0-self.opt['seg_width'])*0.5)
                starty = self.wing_ext_max+self.opt['dist_from_tree']+self.opt['seg_radial_depth']*i

                art2 = None
                lw = v[4] if v[4] else self.opt['seg_line_width']
                if not v[2] or v[2] == 'R':
                    if lw > 0.0:
                        art2 = mpatches.Rectangle(      (startx,starty), width = width, height = height,
                                                        fc = 'none', ec='k', linewidth= lw )
                    art = mpatches.Rectangle(   (startx,starty),
                                                 width = width,
                                                 height = height,
                                                 alpha=v[1],
                                                 color=v[0],
                                                 linewidth= 0.0
                                                 )

                elif v[2] == '^':
                    if lw > 0.0:
                        art2 = mpatches.Polygon( [   [startx,starty],[startx+width/2, starty+height], [startx+width,starty ] ],
                                                     fc = 'none', ec='k', linewidth=  lw)
                    art = mpatches.Polygon( [   [startx,starty],
                                                [startx+width/2, starty+height],
                                                [startx+width,starty ] ],
                                                alpha=v[1], color=v[0],
                                                linewidth=  0.0)
                elif v[2] == 'v':
                    if lw > 0.0:
                        art2 = mpatches.Polygon( [  [startx,starty + height], [startx+width/2, starty], [startx+width,starty + height ] ],
                                                    fc = 'none', ec='k', linewidth=  lw)
                    art = mpatches.Polygon( [   [startx,starty + height],
                                                [startx+width/2, starty],
                                                [startx+width,starty + height ] ],
                                                alpha=v[1], color=v[0],
                                                linewidth=  0.0)
                if art2:
                    ax.add_patch(art2)
                ax.add_patch(art)

                #ax.bar(     clade.pos.rad_min+abs(clade.pos.rad_min-clade.pos.rad_max)*((1.0-self.opt['seg_width'])*0.5),
                #            self.opt['seg_radial_depth'],
                #            width = abs(clade.pos.rad_min-clade.pos.rad_max)*self.opt['seg_width'],
                #            bottom =  self.wing_ext_max+self.opt['dist_from_tree']+self.opt['seg_radial_depth']*i,
                #            alpha=v[1],
                #            color=v[0],
                #            linewidth=  self.opt['seg_line_width'])
        
        def draw_scale():
            nticks = self.opt['branch_length_n_ticks']
            round_prec = 2
            if self.opt['ignore_branch_len']:
                nlev = int(self.__max_br_depth)
                tick_step = 1.0/float(nlev)
                while not ( nlev <= self.opt['branch_length_n_ticks'] 
                            and not int(self.__max_br_depth) % nlev):
                    nlev -= 1
                    tick_step = 1.0/float(nlev)
            else:
                tick_step = (1.0 - self.root.pos.r)/float(nticks)
            tick_pos = np.arange(self.root.pos.r,1.000001,tick_step)
            if self.opt['show_branch_length_ticks']: 
                yticks(tick_pos,[])
            if self.opt['show_branch_length_labels']:
                for t in tick_pos:
                    tv = int(round(t*self.__max_br_depth-1.0,2)) \
                            if self.opt['ignore_branch_len'] \
                            else round(t*self.__max_br_depth,round_prec)
                    if self.tick_labels:
                        tv = self.tick_labels[tv]
                    else: tv = str(tv) 
                    ax.text(    np.pi*0.5,t,
                                tv,
                                ha="center",
                                fontsize = self.opt['branch_length_ticks_font_size'])
        
        def draw_legend():
            ret = []
            h, lleg = ax.get_legend_handles_labels()
            if len(lleg) > 0:
                hl = sorted(zip(h, lleg),key=operator.itemgetter(1))
                h2, l2 = zip(*hl)

                leg = ax.legend(    h2, 
                                    l2, 
                                    bbox_to_anchor=(1.03, 1),
                                    frameon=False, 
                                    prop={'size':self.opt['label_font_size']},
                                    labelspacing=self.opt['legend_label_distance'],
                                    loc=2, 
                                    ncol=self.opt['legend_n_col'],
                                    borderaxespad=0.    )
                ret.append( leg )
            if self.vis_cats:
                labs = []
                for l in self.vis_cats:
                    ll, = ax.bar(   0.0,
                                    0.0,
                                    width=0.0,
                                    bottom=0.0,
                                    alpha=l['alpha'],
                                    color=l['color'],
                                    edgecolor=l['color'],
                                    label=l['label'] )
                    labs.append((ll,l['label']))

                leg2 = ax.legend(   [lp for (lp,la) in labs],
                                    [la for (lp,la) in labs],
                                    bbox_to_anchor=(0.0, 1),
                                    prop={'size':self.opt['category_font_size']},
                                    labelspacing=self.opt['legend_label_distance'],
                                    #                                    prop={'size':self.opt['label_font_size']},
                                    frameon=False,
                                    loc=2,
                                    borderaxespad=0.    )
                ret.append( leg2 )
            if len(lleg) > 0:
                gca().add_artist(leg)
            return ret

        def rec_draw(clade,fcol=self.opt['branch_color']):
            children = clade.get_children()
            if children:
                draw_branch(clade,children,fcol)
            
            rec_col = clade.col if clade.is_highlighted else fcol

            fcls = [rec_draw(c,rec_col) for c in children]

            if not self.opt['draw_taxa']: return clade
            
            for fc in fcls:
                if fc.is_leaf and self.opt['extend_leaves']:
                    draw_extended_leaves(fc)
              
                highlighted_clade = draw_taxa_circle(fc)
                if highlighted_clade:
                    circle_taxa_draw.append(highlighted_clade)

                if fc.is_highlighted:
                    draw_wing(fc)
            
            return clade

        rec_draw(self.root)
        
        
        def rec_draw_sectors(clade,fcol=self.opt['branch_color']):
            children = clade.get_children()
            fcls = [rec_draw_sectors(c) for c in children]
            if not self.opt['draw_taxa']: return clade
            if clade.ext_seg:
                draw_sectors(clade)
            return clade
    

        rec_draw_sectors(self.root)
        coll_b1 = collections.LineCollection(  branch_line_1['frto'],
                                               color = branch_line_1['col'],
                                               linewidths = self.opt['branch_tickness'])
        ax.add_collection(coll_b1)

        coll_b2 = collections.LineCollection(  branch_line_2['frto'],
                                               color = branch_line_2['col'],
                                               linewidths = self.opt['branch_tickness'])
        ax.add_collection(coll_b2)
        
        coll_b3 = collections.LineCollection(   branch_line_3['frto'],
                                                color = branch_line_3['col'],
                                                linewidths = self.opt['branch_tickness'])
        ax.add_collection(coll_b3)
        


        for tc in [taxa_circle, taxa_circle_h]:
            if not tc['rad']: continue
            ax.scatter(     tc['rad'], 
                            tc['r'],
                            #taxa_circle_shape
                            c = tc['col'],
                            s = tc['size'],
                            alpha = 1.0,
                            linewidths = tc['linew'],
                            zorder=12)


        if self.opt['show_branch_length_ticks'] or self.opt['show_branch_length_labels']:
            draw_scale()

        legs = []
        if 'legend_on' not in self.opt or self.opt['legend_on']:
            legs = draw_legend()

        fc = 'w'   
        if 'background_color' in self.opt and self.opt['background_color'] == 'k': 
            def get_col_attr(x):
                return hasattr(x, 'set_color')  and hasattr(x, 'get_color') # and not hasattr(x, 'set_facecolor')

            for o in ax.findobj(get_col_attr):
                col = o.get_color()
                if col == 'k':
                    o.set_color('w')
            fc = 'k'

        a,b = ax.get_ylim()
        ylim((0,b))
        if out_img_file:
            plt.savefig(    out_img_file,
                            dpi=dpi,
                            facecolor=fc,
                            #bbox_inches='tight',
                            #bbox_extra_artists = [l.legendPatch for l in legs],
                            #pad_inches=0.45
                            format = outformat,
                            edgecolor=fc) #,format=self.opt['img_format'])

            plt.close()
        else:
            plt.show()
   
    def read_highlights(self,highlights_file):
        self.labels = {}
        self.label_color = {}
        self.label_cat = {}
        if not highlights_file: 
            return
        with open(highlights_file) as inp_f:
            labels = [l.rstrip().split('\t') 
                    for l in inp_f.readlines() if not l.startswith("#")] 
        for l in labels:
            self.labels[l[0]] = l[1] 
            self.label_cat[l[0]] = l[2]
            if l[3].startswith("_c_"):
                self.label_color[l[0]] = [float(v) for v in l[3].split("_c_")[-1].split("[")[-1].split("]")[0].split(",")]
            else:
                self.label_color[l[0]] = self.colors[l[3]]

    def read_circles(self,circles_file):
        self.cseg = {}
        if not circles_file: 
            return
        with open(circles_file) as inp_f:
            mat = [l.rstrip().split('\t') 
                    for l in inp_f.readlines() if not l.startswith("#")] 
        for m in mat:
            cv = []
            cs = []
            for v in m[1:]:
                v00,bor = v.split("#") if "#" in v else (v,None)
                v0,dep = v00.split("$") if "$" in v00 else (v00,None)
                v1,shape = v0.split("!") if "!" in v0 else (v0,None)
                col,alpha = v1.split(":") if ":" in v1 else [v1,"1.0"]
                
                if col.startswith("_c_"):
                    c = [float(v) for v in col.split("_c_")[-1].split("[")[-1].split("]")[0].split(",")]
                else:
                    c = self.colors[col]
                a = float(alpha)
                cv.append((c,a,shape,float(dep) if dep else None,float(bor) if bor else None))
            self.cseg[m[0]] = cv

    def read_sizes(self,size_file):
        self.sizes = {}
        if not size_file: 
            return
        with open(size_file) as inp_f:
            rows = [l.rstrip().split('\t') 
                    for l in inp_f.readlines() if not l.startswith("#")] 
        for l in rows:
            self.sizes["root."+l[0]] = float(l[1])

    def read_tick_labels(self,ticks_file):
        self.tick_labels = {}
        if not ticks_file:
            return
        with open(ticks_file) as inp_f:
            labels = [l.rstrip().split('\t') 
                        for l in inp_f.readlines() if not l.startswith("#")] 
        for l in labels:
            self.tick_labels[int(l[0])-1] = l[1] 


    default_colors = 'bgrcmy'

    def read_colors(self,colors_file):
        if not colors_file:
            self.opt = {}
            for c in self.default_colors:
                self.opt[c] = c
            return

        self.color_list = []
        self.colors = {}
        with open(colors_file) as inp_f:
            col = [l.rstrip().split('\t') 
                    for l in inp_f.readlines() if not l.startswith("#")]
        for c in col:
            self.color_list.append(c[0])
            self.colors[c[0]] = [float(cc)/255.0 for cc in c[1].split(',')]

    def read_style(self,style_file):
        with open(style_file) as inp_f:
            self.opt = dict([(l.rstrip().split()[0],l.split("#")[0].split()[1:]) 
                                for l in inp_f.readlines() 
                                    if l.strip() and not l.startswith("#")])
        for o in self.opt:
            try:
                v= int(self.opt[o][0])
            except ValueError:
                try:
                    v= float(self.opt[o][0])
                except ValueError:
                    try:
                        v = str(self.opt[o][0])[0]
                    except ValueError:
                        print "not a valid input",self.opt[o][0]
            self.opt[o] = v

