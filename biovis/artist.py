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
        node_type_id =  cell_prop['node_type_id']	
        morphology = morphologies[node_type_id]
        segs_start,segs_end = tr.compute_segs(morphology,cell_prop)

        color = cmap_rgb[cell_prop[color_label]]

        segs_coords, segs_colours = draw_line_segments(segs_start,segs_end,color)
        segments.extend(segs_coords)
        segments_colours.extend(segs_colours)
        
    print "ready to display 3d segments!"

    return segments, segments_colours


def draw_morphologies3D(cells_select_df, morphologies, cmap, color_label, n_faces):
    ''' Draw 3D morphologies'''
    print "...drawing 3D segments..."

    #Test single cell with 3D cylinders
    
    segments3D = []  #collect all segments; convert to numpy array in opengl class; 
    segments3D_colours = []
    segments3D_joints = []
    segments3D_joints_colours = []

    cell_counter = 0
    cmap_rgb = cm.convert_to_rgb(cmap)

    for gid, cell_prop in cells_select_df.iterrows():  
        cell_counter+=1         
        if cell_counter%1000==0: print "... processing cell: ", cell_counter

        #Select colour
        color = cmap_rgb[cell_prop[color_label]]

        #Load segs starts/ends 
        node_type_id =  cell_prop['node_type_id']	
        morphology = morphologies[node_type_id]
        segs_start,segs_end = tr.compute_segs(morphology,cell_prop)

        segs_coords, segs_colours, segs_coords_joints, segs_joints_colours = draw_line_segments3D(segs_start,segs_end,color,n_faces)
        segments3D.extend(segs_coords)
        segments3D_colours.extend(segs_colours)
        segments3D_joints.extend(segs_coords_joints)
        segments3D_joints_colours.extend(segs_joints_colours)

    return segments3D, segments3D_colours, segments3D_joints, segments3D_joints_colours



def draw_line_segments3D(segs_start, segs_end, color, n_faces):

    ''' multiplicate, translate and rotate single segments for plotting of cylinder surfaces 
    '''

    #plotting parameters;
    a1 = 1.; a2 = 1.            #constant radii of cylinders; TODO Load proper dendrite radii
    N = n_faces                 #number of sides in the cylinder; 
    shading_gradient = 0.7      #0.9..0.1, controls the amount of gradient in the shading

    #Make artificial color shading array; it generates a shading gradient as you move around the cylinder faces
    shade_array1 = np.linspace(1.0, shading_gradient, N/2+1); shade_array2 = shade_array1[1:-1][::-1]
    shade_array = np.concatenate((shade_array1, shade_array2), axis=0)

    #convert into an array for multiplication of each segment 
    color_array = [color]*N
    seg_colours = []
    for face, shader in zip(color_array, shade_array):
        seg_colours.append(np.array([face]*4)*shader)
    seg_colours = np.array(seg_colours)
    
    
    #define rotation array
    thetan  = np.arange(N)*2*math.pi/N
    xn=np.cos(thetan)
    zn=np.sin(thetan)
 
    direction_cylinder_axis_prim = [0,1,0]       #primitive direction always up; 
    segs_dl = segs_start - segs_end
    

    #For color normalization by depth; may not be necessary, but looks better with it for now
    #This uses a gradient normalized across y-height of cell; another option is to normalize to depth of column (all cells would thus have same gradient)
    lowest_seg = np.min(segs_start[:,1])
    highest_seg = np.max(segs_start[:, 1])
    range_seg = highest_seg-lowest_seg

    #Required for fusing cylinders together
    segs_array = np.zeros((2*len(segs_dl),3), dtype=np.float32)         #Keep track of segs already parsed; allows proper fusing of segment cylinders 
    nodes_array = np.zeros((2*len(segs_dl),N,3), dtype=np.float32)        #keep track of nodes/cylinders already processed

    segs_coords = []
    segs_joints_coords = []
    segs_colours = []
    segs_joints_colours = []
    for iseg in range(1,segs_dl.shape[0],1):           

        seg_dl = segs_dl[iseg,:]
 
        dlmag = np.linalg.norm(seg_dl)
 
        prim_start=np.zeros((N,3))     
        prim_start[:,0]=a1*xn
        prim_start[:,2]=a1*zn
        
        prim_end=np.zeros((N,3))
        prim_end[:,0]=a2*xn
        prim_end[:,1]=dlmag
        prim_end[:,2]=a2*zn
 
        seg_start = segs_start[iseg,:]
        
        rotv = np.cross(direction_cylinder_axis_prim, seg_dl/dlmag)
        ctheta = np.dot(direction_cylinder_axis_prim, seg_dl/dlmag)

        RotMat = tr.rotation_matrix(rotv, math.pi+math.acos(ctheta))

        nodes_start = seg_start + np.dot(prim_start,RotMat.T)  # use the tranposed matrix since we multiply on the right
        nodes_end = seg_start + np.dot(prim_end,RotMat.T)
        
        #Process to fuse cylinders together; 
        if True: 
            #NB: starts and ends do not always equal each other to float32 precision
            #NB: cylinder faces between sequential segs don't align (unclear why rotations do this); thus, can't blindly fuse them; 
            #Option 1 (implemented): Search for nearby cylinders and use nodes_end from those cylinder ends
                    #problem: cylinders collapse depending on joing angle and can look quite warped.
            #Option 2 : leave cylinders as is but just add additional quad faces betweeen nearby point to close the gaps between cylinders;
                    #NB: store them separately and only merge pre-plotting;
                    #    This is important because synapse plotting relies on seg3d quads being ordered by seg location
            #Option 3: write analytical method for finding midpoints between cylinder face nodes while preserving sizes; will affect synapse placement
            
            ##********* OPTION 1 **********
            #if False: 
                #Fuse_pts = True                         #Fuse cylinders by searching for nearest cylinder face point; works better with higher face count
                #A = (segs_array-segs_start[iseg])       #Compute vector between current seg_start and any existing saved array 
                #B = (A*A).sum(axis=1)**0.5              #Compute euclidean length
                #C = np.argmin(B)                        #Find minimum length, i.e. closest point
                #if B[C]< 1E-3:                          #Determine if min length is below some value close to zero; if so, there is a matching seg nearby to fuse to
                    #if Fuse_pts: 
                        ##Search for each cylinder face point for nearest binding point from nearest cylinder face location; 
                        #for n in range(len(nodes_start)):               #This might be pythonizable; not a priority
                            #A = (nodes_array[C]-nodes_start[n])         #Compute vectors between current seg_start and existing saved array 
                            #B = (A*A).sum(axis=1)**0.5                  #Compute eucledian length
                            #D = np.argmin(B)                            #Find nearest point
                            #nodes_start[n] = nodes_array[C][D]
                    #else:
                        #nodes_start = nodes_array[C]  #Fuse entire cylinder face to previous location; - THIS IS THE CORRECT METHOD IF ROTATIONS OF SEGS WOULD BE ALWAYS ALIGNED
                
                ##Do the same for segs_end
                #A = (segs_array-segs_end[iseg])          
                #B = (A*A).sum(axis=1)**0.5              
                #C = np.argmin(B)                        
                #if B[C]< 1E-3:                          
                    #if Fuse_pts: 
                        #for n in range(len(nodes_end)):
                            #A = (nodes_array[C]-nodes_end[n])       
                            #B = (A*A).sum(axis=1)**0.5             
                            #D = np.argmin(B)   
                            #nodes_end[n] = nodes_array[C][D]
                    #else:
                        #nodes_end = nodes_array[C]

            if True:
                #****** OPTION 2 ***********
                #Do not fuse cylinders, but add additional surfaces to connect nearest nodes on cylinder end faces

                joint_start = []
                joint_end = []       
                
                #Search segs_end for closeby endponints
                A = (segs_array-segs_start[iseg])       #Compute vector between current seg_start and any existing saved array 
                B = (A*A).sum(axis=1)**0.5              #Compute euclidean length
                C = np.argmin(B)                        #Find minimum length, i.e. closest point
                if B[C]< 1E-3:                          #Determine if min length is below some value close to zero; if yes: matching seg nearby

                    #Search for each cylinder face node for nearest point from nearest cylinder face location; 
                    for n in range(len(nodes_start)):               #This might be pythonizable; not a priority
                        A = (nodes_array[C]-nodes_start[n])         #Compute vectors between current seg_start and existing saved array 
                        B = (A*A).sum(axis=1)**0.5                  #Compute eucledian length
                        D = np.argmin(B)                            #Find nearest point
                        joint_start.append(nodes_start[n])
                        joint_end.append(nodes_array[C][D])

                #Do the same for segs_end
                A = (segs_array-segs_end[iseg])          
                B = (A*A).sum(axis=1)**0.5              
                C = np.argmin(B)                        
                if B[C]< 1E-3:   
                    
                    for n in range(len(nodes_end)):               #This might be pythonizable; not a priority
                        A = (nodes_array[C]-nodes_end[n])         #Compute vectors between current seg_start and existing saved array 
                        B = (A*A).sum(axis=1)**0.5                  #Compute eucledian length
                        D = np.argmin(B)                            #Find nearest point
                        joint_end.append(nodes_end[n])
                        joint_start.append(nodes_array[C][D])
                        

        #Keep track of the saved segs and nodes array for use above to fuse/joint creation
        segs_array[iseg*2] = segs_start[iseg];    segs_array[iseg*2+1] = segs_end[iseg]
        nodes_array[iseg*2] = nodes_start;      nodes_array[iseg*2+1] = nodes_end

        #Use nodes to make 3D quads (i.e. faces of cylinders)
        seg_nodes = []                     
        for k in range(N):
            seg_nodes.extend([nodes_start[k], nodes_end[k], nodes_end[(k+1)%N], nodes_start[(k+1)%N]])

        #Use manual shading of segments by depth; ********** SKIP FOR NOW
        #cf = ((nodes_start[0][1]-lowest_seg)/range_seg+2.)/3.       #Normalize colour as a gradient across depth/length of cell
        #color_factor = np.array([cf,cf,cf])
        #seg_colours = seg_colours*color_factor   
    
        #Save nodes/colours to lists
        segs_coords.extend(seg_nodes)
        segs_colours.extend(seg_colours)


        #Same thing, for joints 
        if len(joint_start)>0:
            seg_joints_nodes = []                     
            for k in range(N):
                seg_joints_nodes.extend([joint_start[k], joint_end[k], joint_end[(k+1)%N], joint_start[(k+1)%N]])
            
            segs_joints_coords.extend(seg_joints_nodes)
            segs_joints_colours.extend(seg_colours)            


    
        #if iseg>10: return segs_coords, segs_colours, segs_joints_coords, segs_joints_colours
        
    return segs_coords, segs_colours, segs_joints_coords, segs_joints_colours


def draw_synapses(segments3D, cid, synapses, nodes_df, morphologies, cmap, color_label, n_faces):

    '''segments3D contains the quad surface vertices needed for localization of synapses
    '''

    #Recover the vertices for each segment
    n_segs = len(segments3D)/(n_faces*4)
    print "...n_segs: ", n_segs
    
    segs_vertices = []   
    for p in range(n_segs):
        segs_vertices.append([])
        for k in range(n_faces):
            segs_vertices[p].append([segments3D[p*n_faces*4+k*4], segments3D[p*n_faces*4+k*4+1]])


    sphere_points=[]
    sphere_colours=[]

    cmap_rgb = cm.convert_to_rgb(cmap)          #************** NEED TO HAVE SYNAPSE COLOUR SCHEME AT SOME POINT

    biovis_dir = os.path.dirname(prim.__file__)  # path to package
    sphere_primitive = prim.load_soma_sphere(biovis_dir)

    #Place synapses on dendrites for a particular cell
	    
    node_type_id =  nodes_df.loc[cid]['node_type_id']	
    
#        color = [255,0,0]                         #************* HARDWIRED COLOR SCHEME FOR NOW; REMOVE LATER
    radius = 2                                #************* SYNAPSE SIZE IS FIXED; SHOULD BE SET FROM ANATOMICS

#        for k in range(len(synapses)):
    for seg_id, src_gid in zip(synapses[0],synapses[1]):
        seg_loc = seg_id-1
#        print seg_id
        color_prop = nodes_df.loc[src_gid][color_label]
#            seg_loc = synapses[k]      
#            if seg_loc > n_segs: continue       #************* REMOVE THIS ONCE SYNAPSE FILES MATCH THE CELL AS TO NOT GO OVER MAX #SEGS

        color = cmap_rgb[color_prop]

        #Select a random face and position along face:
        face = np.random.randint(n_faces)
        pos = np.random.rand()
        
        #set 3d location along vertex
        x_synapse = segs_vertices[seg_loc][face][0][0]+ (segs_vertices[seg_loc][face][1][0]- segs_vertices[seg_loc][face][0][0])*pos
        y_synapse = segs_vertices[seg_loc][face][0][1]+ (segs_vertices[seg_loc][face][1][1]- segs_vertices[seg_loc][face][0][1])*pos
        z_synapse = segs_vertices[seg_loc][face][0][2]+ (segs_vertices[seg_loc][face][1][2]- segs_vertices[seg_loc][face][0][2])*pos
                   
        #Compute unit vector along direction away from dendrite cylinder
        #First find vector to opposite side of cylinder same location
        x_synapse_opposite = segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][0]+ (segs_vertices[seg_loc][(face+n_faces/2)%n_faces][1][0]- segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][0])*pos
        y_synapse_opposite = segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][1]+ (segs_vertices[seg_loc][(face+n_faces/2)%n_faces][1][1]- segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][1])*pos
        z_synapse_opposite = segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][2]+ (segs_vertices[seg_loc][(face+n_faces/2)%n_faces][1][2]- segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][2])*pos

        #Second, compute unit vector across cylinder:
        unit_len = np.linalg.norm([x_synapse_opposite-x_synapse, y_synapse_opposite-y_synapse, z_synapse_opposite-z_synapse])
        x_unit = (x_synapse_opposite-x_synapse)/unit_len
        y_unit = (y_synapse_opposite-y_synapse)/unit_len
        z_unit = (z_synapse_opposite-z_synapse)/unit_len
        
        #Third, add vector back in:
        x_synapse -= x_unit*radius
        y_synapse -= y_unit*radius
        z_synapse -= z_unit*radius

        
        syn_centre = np.array([x_synapse,y_synapse,z_synapse])
        
        sphere1_points,sphere1_colours = draw_sphere(radius, syn_centre, color, sphere_primitive)

        sphere_points.extend(sphere1_points)
        sphere_colours.extend(sphere1_colours)


    #TODO: PLACE SYNAPSES ON SOMA
    #Place synapses on soma
    #for gid, cell_prop in cells_select_df.iterrows():  

    print "ready to display synapses"
    
    return sphere_points, sphere_colours
    
    

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
        node_type_id =  cell_prop['node_type_id']    
        morphology = morphologies[node_type_id]
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


def draw_somas(cells_select_df, morphologies, cmap, color_label):

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
    
        node_type_id =  cell_prop['node_type_id']	
        
        radius = 0.5*morphologies[node_type_id]     #Radius comes from first morphology segment   
#        print "...radius: ", radius

        color = cmap_rgb[cell_prop[color_label]]

        sphere1_points,sphere1_colours = draw_sphere(radius,soma_centre,color,sphere_primitive)

        sphere_points.extend(sphere1_points)
        sphere_colours.extend(sphere1_colours)


    print "ready to display spherical somata!"
    
    return sphere_points, sphere_colours



def draw_sphere(radius, soma_centre, color, sphere_primitive):
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



def draw_layers(side,layer_depths, layer_colors, layer_alpha):
    
    layers = []
    color_list = []
    L = side
    for depth,color_name in zip(layer_depths,layer_colors):

        vertex_0 = [-L, -depth, L]
        vertex_1 = [L, -depth, L]
        vertex_2 = [-L, -depth, -L]

        vertex_3 = [-L, -depth, -L]
        vertex_4 = [L, -depth, -L]
        vertex_5 = [L, -depth, L]

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
    
#    frame = -np.array(frame)
#    frame_colours = np.array(frame_colours)

    return frame, frame_colours


def draw_synapses_old(segments3D, cid, syn_df, cells_select_df, morphologies, cmap, color_label, n_faces):

    '''segments3D contains the quad surface vertices needed for localization of synapses
    '''

    #Recover the vertices for each segment
    n_segs = len(segments3D)/(n_faces*4)
    print "...n_segs: ", n_segs
    
    segs_vertices = []   
    for p in range(n_segs):
        segs_vertices.append([])
        for k in range(n_faces):
            segs_vertices[p].append([segments3D[p*n_faces*4+k*4], segments3D[p*n_faces*4+k*4+1]])


    sphere_points=[]
    sphere_colours=[]

    cmap_rgb = cm.convert_to_rgb(cmap)          #************** NEED TO HAVE SYNAPSE COLOUR SCHEME AT SOME POINT

    biovis_dir = os.path.dirname(prim.__file__)  # path to package
    sphere_primitive = prim.load_soma_sphere(biovis_dir)

    #Place synapses on dendrites
    for gid, cell_prop in cells_select_df.iterrows():  
    
        node_type_id =  cell_prop['node_type_id']	
        
#        color = [255,0,0]                         #************* HARDWIRED COLOR SCHEME FOR NOW; REMOVE LATER
        radius = 1                                #************* SYNAPSE SIZE IS FIXED; SHOULD BE SET FROM ANATOMICS

#        for k in range(len(synapses)):
        for syn_id, syn_prop in syn_df.iterrows():
            seg_loc = syn_prop["seg_id"]
            color_label = syn_prop["src_label"]
#            seg_loc = synapses[k]      
#            if seg_loc > n_segs: continue       #************* REMOVE THIS ONCE SYNAPSE FILES MATCH THE CELL AS TO NOT GO OVER MAX #SEGS

            color = cmap_rgb[color_label]

            #Select a random face and position along face:
            face = np.random.randint(n_faces)
            pos = np.random.rand()
            
            #set 3d location along vertex
            x_synapse = segs_vertices[seg_loc][face][0][0]+ (segs_vertices[seg_loc][face][1][0]- segs_vertices[seg_loc][face][0][0])*pos
            y_synapse = segs_vertices[seg_loc][face][0][1]+ (segs_vertices[seg_loc][face][1][1]- segs_vertices[seg_loc][face][0][1])*pos
            z_synapse = segs_vertices[seg_loc][face][0][2]+ (segs_vertices[seg_loc][face][1][2]- segs_vertices[seg_loc][face][0][2])*pos
                       
            #Compute unit vector along direction away from dendrite cylinder
            #First find vector to opposite side of cylinder same location
            x_synapse_opposite = segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][0]+ (segs_vertices[seg_loc][(face+n_faces/2)%n_faces][1][0]- segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][0])*pos
            y_synapse_opposite = segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][1]+ (segs_vertices[seg_loc][(face+n_faces/2)%n_faces][1][1]- segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][1])*pos
            z_synapse_opposite = segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][2]+ (segs_vertices[seg_loc][(face+n_faces/2)%n_faces][1][2]- segs_vertices[seg_loc][(face+n_faces/2)%n_faces][0][2])*pos

            #Second, compute unit vector across cylinder:
            unit_len = np.linalg.norm([x_synapse_opposite-x_synapse, y_synapse_opposite-y_synapse, z_synapse_opposite-z_synapse])
            x_unit = (x_synapse_opposite-x_synapse)/unit_len
            y_unit = (y_synapse_opposite-y_synapse)/unit_len
            z_unit = (z_synapse_opposite-z_synapse)/unit_len
            
            #Third, add vector back in:
            x_synapse -= x_unit*radius
            y_synapse -= y_unit*radius
            z_synapse -= z_unit*radius

            
            syn_centre = np.array([x_synapse,y_synapse,z_synapse])
            
            sphere1_points,sphere1_colours = draw_sphere(radius, syn_centre, color, sphere_primitive)

            sphere_points.extend(sphere1_points)
            sphere_colours.extend(sphere1_colours)


    #TODO: PLACE SYNAPSES ON SOMA
    #Place synapses on soma
    #for gid, cell_prop in cells_select_df.iterrows():  

    print "ready to display synapses"
    
    return sphere_points, sphere_colours


