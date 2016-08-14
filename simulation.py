#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import os
import glob
import numpy as np
#import struct
#import string, re
#import scipy
#import cPickle as pickle
#import gc
#from skimage.measure import block_reduce
#import shutil

import vismorph as vm 
import pandas as pd
import pylab as plt


class Simulation(object):      
    
    def __init__(self, name, home_dir):
        
        self.name = name
        self.home_dir = home_dir

        self.load_nodes()
        #self.load
    
    def load_nodes(self):
        print "...loading cells...",
    
        node_file_name = self.home_dir+self.name + "/nodes/cells.csv"
        node_prop_file_name = self.home_dir + self.name + "/nodes/cm_perisomatic.csv"
        morph_dir = self.home_dir+self.name + "/morph_segs"

        c_df = pd.read_csv(node_file_name, sep=' ')
        c_df.set_index('id',inplace=True)

        cm_df = pd.read_csv(node_prop_file_name, sep=' ')
        cm_df.set_index('model_id',inplace=True)

        ncells = len(c_df.index) # total number of simulated cells
        print "...total cells: ", ncells

        self.cells_prop_df = pd.merge(left=c_df,
                                right=cm_df, 
                                how='left', 
                                left_on='model_id', 
                                right_index=True) # use 'model_id' key to merge, for right table the "model_id" is an index


        self.morphologies = vm.load_morphologies(cm_df, morph_dir)
        print "...# morphologies: ", len(self.morphologies)          #*****************Cat: not sure what this is; is it # of unique morphologies? 


    def load_query(self, main_widget):
        print "...selecting cells to display..."
        
        selection = main_widget.query_string.text() 
        # query examples:
        #selection = "cre_line=='Rorb' & location=='VisL4' & (z_soma**2+x_soma**2)<200**2"
        #selection = "location=='VisL4' & (z_soma**2+x_soma**2)<300**2"
        #selection = "ei=='excitatory'"
        #selection = "model_id==314642645"

        cells_select_df = self.cells_prop_df.query(selection)
        ncells_select = len(cells_select_df.index)
        print "number of cells satisfying selection:", ncells_select

        #nskip=200 # skip cells when there are too many to display
        nskip = float(main_widget.rarifier.text())
        if nskip>100: print "... can't display more than 100%..."; return
        rows_display = range(0, ncells_select, int(100./nskip))
        if len(rows_display)<1: print "... zero cells to display..."; return
        
        self.cells_display_df = cells_select_df.iloc[rows_display]
        print "number of cells to display:", len(self.cells_display_df.index)
        
        
    def show_query(self):
        
        vm.cell_morphs3d(self.cells_display_df, self.morphologies)
        
        #vm.cell_somata3d(self.cells_display_df)

        plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
