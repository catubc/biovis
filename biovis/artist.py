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

def draw_morphologies(cells_select_df, morphologies, cmap,color_label):
    

    segments = []  #collect all segments; convert to numpy array in opengl class; 
    segments_colours = []

    print "... processing cell: ",
    cell_counter = 0
    cmap_rgb = cm.convert_to_rgb(cmap)

    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%1000==0: print cell_counter,
        model_id =  cell_prop['model_id']	
        morphology = morphologies[model_id]
        segs_start,segs_end = tr.compute_segs(morphology,cell_prop)

        color = cmap_rgb[cell_prop[color_label]]

        segs_coords,segs_colours = draw_line_segments(segs_start,segs_end,color)
        segments.extend(segs_coords)
        segments_colours.extend(segs_colours)
        
    print "ready to display 3d segments!"

    return segments, segments_colours


def draw_morphologies3D(cells_select_df, morphologies, cmap, color_label):
    ''' Draw 3D morphologies'''
    print "...drawing 3D segments..."

    #Test single cell with 3D cylinders
    
    segments3D = []  #collect all segments; convert to numpy array in opengl class; 
    segments3D_colours = []

    cell_counter = 0
    cmap_rgb = cm.convert_to_rgb(cmap)

    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%1000==0: print "... processing cell: ", cell_counter

        #Select colour
        color = cmap_rgb[cell_prop[color_label]]

        #Load segs starts/ends 
        model_id =  cell_prop['model_id']	
        morphology = morphologies[model_id]
        segs_start,segs_end = tr.compute_segs(morphology,cell_prop)

        #Make line seg arrays for opengl
        segs_coords, segs_colours = draw_line_segments3D(segs_start,segs_end,color)
        segments3D.extend(segs_coords)
        segments3D_colours.extend(segs_colours)
        
        break       #************************************DEBUG ONLY; Process single cell  


    return segments3D, segments3D_colours



def draw_line_segments3D(starts, ends, color):

    '''multiplicate, translate and rotate single segments for plotting of cylinder surfaces 
    '''

    print "... # of segs: ", len(starts)
    #print starts
    #initialize holding matrices
    segs_start = []
    segs_end = []

    #make 6 sided cylinders
    a1 = 10.; a2 = 10.      #constant radii of cylinders
    theta_n = []
    for k in range(6):  
        theta_n.append(2*math.pi/6.*k)


    #loop over all segments and tranform the primitive along correct angles
    #for l in range(len(starts)): 
    for l in range(2):  #Do just a few segments while debugging
        h = np.linalg.norm(ends[l]-starts[l])
        print "...length of seg: ", l, " = ", h

        #Save locations of primitve starts/ends to be rotated below
        temp_starts = []
        temp_ends = []
        for k in range(len(theta_n)):       #Is this correct?!
            temp_starts.append([a1*math.cos(theta_n[k]), 0, a1*math.sin(theta_n[k])])
            temp_ends.append([a2*math.cos(theta_n[k]), h, a2*math.sin(theta_n[k])])

        #rotate 6 lines
        for k in range(len(temp_starts)):
            direction_cylinder_axis = [0,1,0]       #primitive direciton always up; is this correct ?!
            direction_segment_center = [temp_ends[k][0]-temp_starts[k][0],      #Swapping y-axis <-> z-axis; correct?! 
                                        temp_ends[k][2]-temp_starts[k][2], 
                                        temp_ends[k][1]-temp_starts[k][1]]
                        
            #************getting lost around this place****************
            rotv = np.cross(direction_cylinder_axis, direction_segment_center)
            ctheta = np.dot(direction_cylinder_axis, direction_segment_center)

            r = tr.rotation_matrix(rotv, math.acos(ctheta))
            
            segs_start.append(np.dot(temp_starts[k], r)+starts[l])
            segs_end.append(np.dot(temp_ends[k], r)+ends[l])


    #make coords arrays
    segs_coords = np.empty((len(segs_start)*2, 3), dtype=np.float32)
    segs_coords[::2] = segs_start; segs_coords[1::2] = segs_end

    segs_colours = [color]*len(segs_coords)

    return segs_coords, segs_colours


def draw_line_segments(segs_start,segs_end,color):
    '''
    interleave starts/ends into arrays for plotting;
    '''
    segs_coords = np.empty((len(segs_start)*2, 3), dtype=np.float32)
    segs_coords[::2] = segs_start; segs_coords[1::2] = segs_end

    segs_colours = [color]*len(segs_coords)

    return segs_coords, segs_colours



def draw_slice(cells_select_df, morphologies, cmap,color_label,xplane_range):
    

    segments = []  #collect all segments; convert to numpy afterwards
    segments_colours = []

    print "... processing cell: ",
    cell_counter = 0
    cmap_rgb = cm.convert_to_rgb(cmap)

    xmin = xplane_range[0]
    xmax = xplane_range[1]
    
    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%1000==0: print cell_counter,
        model_id =  cell_prop['model_id']    
        morphology = morphologies[model_id]
        segs_start,segs_end = tr.compute_segs(morphology,cell_prop)

        segs_center = 0.5*(segs_start+segs_end)
        
        ix_slice_max = np.where(segs_center[:,0]<xmax)
        ix_slice_min = np.where(segs_center[:,0]>xmin)
        ix_slice = np.intersect1d(ix_slice_min,ix_slice_max)
        
        color = cmap_rgb[cell_prop[color_label]]

        segs_start_slice = np.squeeze(segs_start[ix_slice,:])
        segs_end_slice = np.squeeze(segs_end[ix_slice,:])
        
        segs_coords,segs_colours = draw_line_segments(segs_start_slice,segs_end_slice,color)
        segments.extend(segs_coords)
        segments_colours.extend(segs_colours)
        
    print "ready to display 3d segments!"

    return segments, segments_colours


def draw_somas(cells_select_df, soma_sizes, cmap, color_label):


    sphere_points=[]
    sphere_colours=[]

    cell_counter = 0
    cmap_rgb = cm.convert_to_rgb(cmap)

    biovis_dir = os.path.dirname(prim.__file__)  # path to package
    sphere_primitive = prim.load_soma_sphere(biovis_dir)

    print "... processing cell: ",
    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%1000==0: print cell_counter,
        
        x_soma = cell_prop['x_soma'] # needed to do those in sequence because Series did not inherite the dtype from DataFrame 
        y_soma = cell_prop['y_soma'] # need to look for a more elegant solution to this 
        z_soma = cell_prop['z_soma']
        soma_centre = np.array([x_soma,y_soma,z_soma])
    
        model_id =  cell_prop['model_id']	
        
        radius = 0.5*soma_sizes[model_id]       
        color = cmap_rgb[cell_prop[color_label]]

        sphere1_points,sphere1_colours = draw_sphere(radius,soma_centre,color,sphere_primitive)

        sphere_points.extend(sphere1_points)
        sphere_colours.extend(sphere1_colours)


    print "ready to display spherical somata!"
    
    return sphere_points, sphere_colours



def draw_sphere(radius, soma_centre,color,sphere_primitive):
    '''
        #SLG: please optimize
    '''
    sphere1_points = []
    sphere1_colours =[]


    vertices, triangle_faces, lowest_vertex_triangle = sphere_primitive

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



def draw_layers(layer_depths, layer_colors, layer_alpha):
    
    layers = []
    color_list = []
    for depth,color_name in zip(layer_depths,layer_colors):

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
        
        col_rgb = cm.colors_lib[color_name]
        color_list.extend([col_rgb]*6)    
     
    layers = np.float32(layers)
    color_list = np.float32(color_list)
    
    layers_colours = np.insert(color_list, 3, layer_alpha, axis=1)
      
    return layers, layers_colours



def draw_frame(box_coords, box_colour):
    
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

