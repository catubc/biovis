#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License


import h5py
import pandas as pd
import numpy as np

#********************************************************************

def load_nodes(cells_file_name, cell_models_file_name):
    print "...importing nodes: ",


    c_df = pd.read_csv(cells_file_name, sep=' ')
    c_df.set_index('id',inplace=True)

    cm_df = pd.read_csv(cell_models_file_name, sep=' ')
    cm_df.set_index('model_id',inplace=True)

    ncells = len(c_df.index) # total number of simulated cells

    nodes_df = pd.merge(left=c_df,
                            right=cm_df, 
                            how='left', 
                            left_on='model_id', 
                            right_index=True) # use 'model_id' key to merge, for right table the "model_id" is an index
    
    print ncells
    return nodes_df, cm_df

   
    
def load_morphologies(morph_dir, cm_df):

    morphologies = {}
    soma_sizes = {}

    for model_id, morph_prop in cm_df.iterrows():
    
        file_name = morph_dir+'/%d.h5' % (model_id)
    
        f5 = h5py.File(file_name,'r')
        segs_start = f5['segs_start'][...]
        segs_end = f5['segs_end'][...]
    
        morphologies[model_id] = {}
    
        morphologies[model_id]["segs_start"] = segs_start
        morphologies[model_id]["segs_end"] = segs_end

        soma_start = segs_start[0]
        soma_end =   segs_end[0]

        size=np.linalg.norm(soma_start-soma_end)          #Size of cell soma;
        soma_sizes[model_id] = size
        

        f5.close()

    print "# unique morphologies: ", len(morphologies)

        
    return morphologies,soma_sizes



def load_synapses(synapses_dir, gid):
    
    import csv
    with open(synapses_dir+'/'+str(gid)+'.csv', 'rb') as csvfile:
        rows = csv.reader(csvfile)
        next(rows, None)    #Skip header
        syn_locs = []
        for row in rows:
            syn_locs.append(row[0].split()[1])
    syn_locs = np.int16(syn_locs)
    
    return syn_locs
    
