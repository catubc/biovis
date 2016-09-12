#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import os
import numpy as np
import h5py
import math
import csv

import pandas as pd

#****************** COLOUR DICTIONARY ****************************

RED = 255, 0, 0
ORANGE = 255, 127, 0
YELLOW = 255, 255, 0
GREEN = 0, 255, 0
CYAN = 0, 255, 255
LIGHTBLUE = 0, 127, 255
BLUE = 0, 0, 255
VIOLET = 127, 0, 255
MAGENTA = 255, 0, 255  
GREY = 85, 85, 85
WHITE = 255, 255, 255
DARK_GREY = 30, 30, 30
PURPLE = 154, 44, 209
CMAP = np.array([RED, ORANGE, YELLOW, GREEN, CYAN, LIGHTBLUE, BLUE, VIOLET, MAGENTA,
                 GREY, WHITE, PURPLE], dtype=np.uint8)


#********************************************************************

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
    


def plot_morphologies(cells_select_df, morphologies):
    

    segments = []  #collect all segments; convert to numpy afterwards
    for gid, cell_prop in cells_select_df.iterrows():  
        if gid%1000==0: print "... processing cell: ", gid
        
        x_soma = cell_prop['x_soma'] # needed to do those in sequence because Series did not inherite the dtype from DataFrame 
        y_soma = cell_prop['y_soma'] # need to look for a more elegant solution to this 
        z_soma = cell_prop['z_soma']
        pos_soma = np.array([x_soma,y_soma,z_soma])
    
        phi_y = cell_prop['rotation_angle_yaxis']
        phi_z = cell_prop['rotation_angle_zaxis']
    
        RotY = rotation_matrix([0,1,0],phi_y)
        RotZ = rotation_matrix([0,0,1],-phi_z)
        RotYZ = RotY.dot(RotZ)	# apply two rotations

        model_id =  cell_prop['model_id']	
        morph_segs_start = morphologies[model_id]["segs_start"]
        morph_segs_end =   morphologies[model_id]["segs_end"] 
            
        segs_start = pos_soma + np.dot(morph_segs_start,RotYZ.T)  # use the tranposed matrix since we multiply on the right
        segs_end = pos_soma + np.dot(morph_segs_end,RotYZ.T)
       
        tot_segs = np.empty((len(segs_start)*2, 3), dtype=np.float32)
        tot_segs[::2]=segs_start; tot_segs[1::2]=segs_end
                    
        segments.extend(tot_segs)

    print "... done loading cells..."
         
    print "...generating cell colours..."
    segments_colours = []
    for k in range(len(segments)):     #should do this in one step...
        segments_colours.append([0,255,0])
    
    print "...converting lists to arrays ..."
    segments = np.array(segments)
    segments_colours = np.array(segments_colours)

    print "... done all ..."

    return segments, segments_colours



def plot_somas(cells_select_df, morphologies):

    vertices, triangle_faces, lowest_vertex_triangle = load_soma_sphere()    #Open sphere primitive

    sphere_points=[]
    sphere_colours=[]

    for gid, cell_prop in cells_select_df.iterrows():  
        if gid%1000==0: print "... processing cell: ", gid
        
        x_soma = cell_prop['x_soma'] # needed to do those in sequence because Series did not inherite the dtype from DataFrame 
        y_soma = cell_prop['y_soma'] # need to look for a more elegant solution to this 
        z_soma = cell_prop['z_soma']
        pos_soma = np.array([x_soma,y_soma,z_soma])
    
        phi_y = cell_prop['rotation_angle_yaxis']
        phi_z = cell_prop['rotation_angle_zaxis']
    
        RotY = rotation_matrix([0,1,0],phi_y)
        RotZ = rotation_matrix([0,0,1],-phi_z)
        RotYZ = RotY.dot(RotZ)	# apply two rotations
            
        model_id =  cell_prop['model_id']	
        
        soma_start = morphologies[model_id]["segs_start"][0]
        soma_end =   morphologies[model_id]["segs_end"][0]
        
        soma_start = pos_soma + np.dot(soma_start,RotYZ.T)
        soma_end = pos_soma + np.dot(soma_end,RotYZ.T)
        
       
        size=np.linalg.norm(soma_start-soma_end)          #Size of cell soma;
        
        soma_centre=(soma_end+soma_start)/2.


        #Make triangle surfaces
        for j in range(len(triangle_faces)):
            sphere_points.append([[vertices[int(triangle_faces[j][0])-1][0]*size/2.+soma_centre[0],
                                        vertices[int(triangle_faces[j][0])-1][1]*size/2.+soma_centre[1]-size/2,
                                        vertices[int(triangle_faces[j][0])-1][2]*size/2.+soma_centre[2]],
                                        
                                         [vertices[int(triangle_faces[j][1])-1][0]*size/2.+soma_centre[0], 
                                        vertices[int(triangle_faces[j][1])-1][1]*size/2.+soma_centre[1]-size/2, 
                                        vertices[int(triangle_faces[j][1])-1][2]*size/2.+soma_centre[2]], 
                                        
                                         [vertices[int(triangle_faces[j][2])-1][0]*size/2.+soma_centre[0],
                                         vertices[int(triangle_faces[j][2])-1][1]*size/2.+soma_centre[1]-size/2,
                                         vertices[int(triangle_faces[j][2])-1][2]*size/2.+soma_centre[2]]                                             
                                         ])

                                
            color_factor = ((lowest_vertex_triangle[j])+2.)/4. #Cat: take z direction and normalize
            #self.triangle_colours.append([[int(CMAP[gid%10][0]*color_factor),
            #                                int(CMAP[gid%10][1]*color_factor),
            #                                int(CMAP[gid%10][2]*color_factor)
            #                                ]]*3)
             
            for k in range(3):
                sphere_colours.append([int(CMAP[gid%10][0]*color_factor),
                                              int(CMAP[gid%10][1]*color_factor),
                                              int(CMAP[gid%10][2]*color_factor)
                                              ])
    
    
    sphere_points = np.array(sphere_points)
    sphere_colours = np.array(sphere_colours)

    return sphere_points, sphere_colours

def load_soma_sphere():     

    #Load sphere vertices
    f = open('static/sphere_vertices.csv', 'rt')
    vertices_list = list(csv.reader(f))
    f.close()
    vertices = []
    for row in vertices_list:
        vertices.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3])])
    
    #Load sphere faces
    f = open('static/sphere_faces.csv', 'rt')
    face_list = list(csv.reader(f))
    f.close()
    triangle_faces=[]
    quad_faces=[]
    for row in face_list:
        if (len(row[0].split(" ")[1:]))==4:     #Can use quad spheres or triangle spheres
            quad_faces.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3]), float(row[0].split(" ")[4])])
        else:
            triangle_faces.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3])])


    lowest_vertex_triangle = []
    for j in range(len(triangle_faces)):
        lowest_vertex_triangle.append(min(vertices[int(triangle_faces[j][0])-1][1],
                            vertices[int(triangle_faces[j][1])-1][1],
                            vertices[int(triangle_faces[j][2])-1][1]))

    #lowest_vertex_quad = []
    #for j in range(len(self.quad_faces)):
        #lowest_vertex_quad.append(min(self.vertices[int(self.quad_faces[j][0])-1][1],
                #self.vertices[int(self.quad_faces[j][1])-1][1],
                #self.vertices[int(self.quad_faces[j][2])-1][1],
                #self.vertices[int(self.quad_faces[j][3])-1][1]))
    #self.lowest_vertex_quad=lowest_vertex_quad


    return vertices, triangle_faces, lowest_vertex_triangle

def rotation_matrix( axis, theta):
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



def load_layers(layer_depths, layer_colors, layer_alpha):
    
    layers = []

    for depth in layer_depths:
        vertex_0 = [-500., -depth, 500.]
        vertex_1 = [500., -depth, 500.]
        vertex_2 = [-500., -depth, -500.]

        vertex_3 = [-500., -depth, -500.]
        vertex_4 = [500., -depth, -500.]
        vertex_5 = [500., -depth, 500.]

        layers.append(vertex_0)
        layers.append(vertex_1)
        layers.append(vertex_2)
        layers.append(vertex_3)
        layers.append(vertex_4)
        layers.append(vertex_5)
            
    layers_colours = layer_colors * len(layer_depths)*6  #uint8; Need colour for every node, not every vertex; i.e. 2 x no. vertices  

    layers = np.float32(layers)
    layers_colours = np.float32(layers_colours)
    
    layers_colours = np.insert(layers_colours, 3, layer_alpha, axis=1)
      
    return layers, layers_colours


def load_frame(box_coords, box_colour):
    
    frame = []
    frame.append(box_coords[0])
    frame.append(box_coords[1])
    frame.append(box_coords[1])
    frame.append(box_coords[2])
    frame.append(box_coords[2])
    frame.append(box_coords[3])
    frame.append(box_coords[3])
    frame.append(box_coords[0])

    frame.append(box_coords[4])
    frame.append(box_coords[5])
    frame.append(box_coords[5])
    frame.append(box_coords[6])
    frame.append(box_coords[6])
    frame.append(box_coords[7])
    frame.append(box_coords[7])
    frame.append(box_coords[4])

    frame.append(box_coords[4])
    frame.append(box_coords[0])
    frame.append(box_coords[5])
    frame.append(box_coords[1])
    frame.append(box_coords[6])
    frame.append(box_coords[2])
    frame.append(box_coords[7])
    frame.append(box_coords[3])

    frame_colours = box_colour*72  
    
    frame = -np.array(frame)
    frame_colours = np.array(frame_colours)

    return frame, frame_colours

