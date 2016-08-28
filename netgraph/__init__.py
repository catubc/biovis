""" initiatialization file
"""

import os
import glob
import numpy as np
import h5py
import math

import pandas as pd
import pylab as plt


def load_nodes(cells_file_name, cell_models_file_name):
    print "...importing nodes..."

    #node_file_name = self.home_dir+self.name + "/nodes/cells.csv"
    #node_prop_file_name = self.home_dir + self.name + "/nodes/cm_perisomatic.csv"
    #morph_dir = self.home_dir+self.name + "/morph_segs"

    c_df = pd.read_csv(cells_file_name, sep=' ')
    c_df.set_index('id',inplace=True)

    cm_df = pd.read_csv(cell_models_file_name, sep=' ')
    cm_df.set_index('model_id',inplace=True)

    ncells = len(c_df.index) # total number of simulated cells
    print "...total # cells simulated: ", ncells

    cells_prop_df = pd.merge(left=c_df,
                            right=cm_df, 
                            how='left', 
                            left_on='model_id', 
                            right_index=True) # use 'model_id' key to merge, for right table the "model_id" is an index
    
    return cm_df, cells_prop_df

#def load_morph_segs(file_name_morph)
#
#    morphologies = load_morphologies(cm_df, )
#    print "...# morphologies: ", len(morphologies)          #*****************Cat: not sure what this is; is it # of unique morphologies? 

   
    
def load_morphologies(morph_dir, cm_df):

    morphologies = {}
    
    for model_id, morph_prop in cm_df.iterrows():
    
        file_name = morph_dir+'/%d.h5' % (model_id)
    
        f5 = h5py.File(file_name,'r')
        segs_start = f5['segs_start'][...]
        segs_end = f5['segs_end'][...]
    
        morphologies[model_id] = {}
    
        morphologies[model_id]["segs_start"] = segs_start
        morphologies[model_id]["segs_end"] = segs_end
        f5.close()
        
    return morphologies
    



def set_axes_equal(self, ax):
    '''
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]; x_mean = np.mean(x_limits)
    y_range = y_limits[1] - y_limits[0]; y_mean = np.mean(y_limits)
    z_range = z_limits[1] - z_limits[0]; z_mean = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])
    
    

def cell_morphs3d(self):

    fig = plt.figure(232)
    ax = fig.gca(projection='3d')
    
    col='k'
    
    RotX = np.array([[1, 0, 0],    # rotate around x axis
                     [0, 0, 1],
                     [0, 1, 0]])
    
    
    for gid,cell_prop in self.cells_display_df.iterrows():  
    
    
        x_soma = cell_prop['x_soma'] # needed to do those in sequence because Series did not inherite the dtype from DataFrame 
        y_soma = cell_prop['y_soma'] # need to look for a more elegant solution to this 
        z_soma = cell_prop['z_soma']
        pos_soma = np.array([x_soma,y_soma,z_soma])
    
        phi_y = cell_prop['rotation_angle_yaxis']
        phi_z = cell_prop['rotation_angle_zaxis']
    
        RotY = self.rotation_matrix([0,1,0],phi_y)
        RotZ = self.rotation_matrix([0,0,1],-phi_z)
        RotYZ = RotY.dot(RotZ)	# apply two rotations
    
    
        model_id =  cell_prop['model_id']	
        morph_segs_start = self.morphologies[model_id]["segs_start"]
        morph_segs_end =   self.morphologies[model_id]["segs_end"] 
    
    
        segs_start = pos_soma + np.dot(morph_segs_start,RotYZ.T)  # use the tranposed matrix since we multiply on the right
        segs_end = pos_soma + np.dot(morph_segs_end,RotYZ.T)
        pos_soma = np.dot(pos_soma,RotX.T)*1E-3
    
    
        segs_start = np.dot(segs_start,RotX.T)*1E-3	
        segs_end = np.dot(segs_end,RotX.T)*1E-3		# for display purposes switch z and y, convert to mm:
    
        col_segs = self.build_collection(segs_start,segs_end)
    
        ax.add_collection(col_segs)	# plot segments
    
        ax.scatter(pos_soma[0], pos_soma[1], pos_soma[2],color='red',alpha=0.4) # plot somata
    
    
    self.set_axes_equal(ax)
    ax.set_aspect('equal')
    
    ax.set_xlim([-0.4,0.4])
    ax.set_zlim([-1,0.2])
    ax.set_ylim([-0.4,0.4])
    
    
    plt.locator_params(nbins=4)
    
    ax.set_xlabel('lat')
    ax.set_ylabel('lat')
    ax.set_zlabel('depth')

    
def cell_somata3d(self):

    fig = plt.figure(233)
    ax = fig.gca(projection='3d')
    
    color='r'
    marker = 'o'
    
    x = self.cell_display_df['x_soma']*1E-3
    y = self.cell_display_df['y_soma']*1E-3
    z = self.cell_display_df['z_soma']*1E-3
    
    ax.scatter(x,z,y,color=color, alpha=0.4,marker=marker)
    
    self.set_axes_equal(ax)
    ax.set_aspect('equal')
    
    ax.set_xlim([-0.4,0.4])
    ax.set_zlim([-1,0.2])
    ax.set_ylim([-0.4,0.4])
    
    plt.locator_params(nbins=4)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_zlabel('y')

    
    
def rotation_matrix(self, axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d

    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        
    

def build_collection(self, segs_start,segs_end):

    mysegs = []

    for seg_start,seg_end in zip(segs_start,segs_end):
        mysegs.append((seg_start,seg_end))

    col_segs = mplot3d.art3d.Line3DCollection(mysegs,colors='k',alpha=0.2)

    return col_segs