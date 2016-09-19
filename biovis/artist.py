#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import os
import numpy as np
import csv

import numpy as np
import math
import primitives as prim
import transformer as tr

import cmaps as cm

def construct_morphologies(cells_select_df, morphologies, cmap,color_label):
    

    segments = []  #collect all segments; convert to numpy afterwards
    print "... processing cell: ",
    cell_counter = 0
    segments_colours = []
    cmap_rgb = cm.convert_to_rgb(cmap)

    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%100==0: print cell_counter,
        model_id =  cell_prop['model_id']	
        morphology = morphologies[model_id]
        segs_start,segs_end = tr.get_segs(morphology,cell_prop)

        tot_segs = np.empty((len(segs_start)*2, 3), dtype=np.float32)
        tot_segs[::2]=segs_start; tot_segs[1::2]=segs_end
        
        segments.extend(tot_segs)
        color = cmap_rgb[cell_prop[color_label]]
        cell_seg_colors = [color]*len(tot_segs)
        segments_colours.extend(cell_seg_colors)
        
    print "ready to display 3d segments!"

    return segments, segments_colours



def construct_somas(cells_select_df, soma_sizes, cmap, color_label):

    vertices, triangle_faces, lowest_vertex_triangle = prim.load_soma_sphere()    #Open sphere primitive

    sphere_points=[]
    sphere_colours=[]

    cell_counter = 0
    cmap_rgb = cm.convert_to_rgb(cmap)

    print "... processing cell: ",
    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%100==0: print cell_counter,
        
        x_soma = cell_prop['x_soma'] # needed to do those in sequence because Series did not inherite the dtype from DataFrame 
        y_soma = cell_prop['y_soma'] # need to look for a more elegant solution to this 
        z_soma = cell_prop['z_soma']
        soma_centre = np.array([x_soma,y_soma,z_soma])
    
        model_id =  cell_prop['model_id']	
        
        radius = 0.5*soma_sizes[model_id]       
        color = cmap_rgb[cell_prop[color_label]]


        sphere1_points,sphere1_colours = construct_sphere(radius,soma_centre,color)

        sphere_points.extend(sphere1_points)
        sphere_colours.extend(sphere1_colours)


    print "ready to display!"
    
    return sphere_points, sphere_colours



def construct_sphere(radius, soma_centre,color):
    '''
        #SLG: please optimize
    '''
    sphere1_points = []
    sphere1_colours =[]

    vertices = prim.sphere['vertice']
    triangle_faces = prim.sphere['triangle_faces']
    lowest_vertex_triangle = prim.sphere['lowest_vertex_triangle']


    for j in range(len(triangle_faces)):
        sphere1_points.append([[vertices[int(triangle_faces[j][0])-1][0]*radius+soma_centre[0],
                                    vertices[int(triangle_faces[j][0])-1][1]*radius+soma_centre[1]-radius,
                                    vertices[int(triangle_faces[j][0])-1][2]*radius+soma_centre[2]],
                                       
                                     [vertices[int(triangle_faces[j][1])-1][0]*radius+soma_centre[0], 
                                    vertices[int(triangle_faces[j][1])-1][1]*radius+soma_centre[1]-radius, 
                                    vertices[int(triangle_faces[j][1])-1][2]*radius+soma_centre[2]], 
                                       
                                     [vertices[int(triangle_faces[j][2])-1][0]*radius+soma_centre[0],
                                     vertices[int(triangle_faces[j][2])-1][1]*radius+soma_centre[1]-radius,
                                     vertices[int(triangle_faces[j][2])-1][2]*radius+soma_centre[2]]                                             
                                     ])

        color_factor = ((lowest_vertex_triangle[j])+2.)/4. #Cat: take z direction and normalize

        for k in range(3):
            sphere1_colours.append([int(color[0]*color_factor),
                                        int(color[1]*color_factor),
                                        int(color[2]*color_factor)
                                        ])


    return sphere1_points,sphere1_colours


def construct_layers(layer_depths, layer_colors, layer_alpha):
    
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


def construct_frame(box_coords, box_colour):
    
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

