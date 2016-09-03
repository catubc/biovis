import numpy as np

def rotation_matrix(axis, theta):
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



def map_dendrites(self, cells_select_df):
    
    #RotX = np.array([[1, 0, 0],    # rotate around x axis      #Cat: NOT USED
    #             [0, 0, 1],
    #             [0, 1, 0]])

    self.segments = []  #collect all segments; convert to numpy afterwards
    for gid, cell_prop in cells_select_df.iterrows():  
        if gid%1000==0: print "... processing cell: ", gid
        
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
       
        tot_segs = np.empty((len(segs_start)*2, 3), dtype=np.float32)
        tot_segs[::2]=segs_start; tot_segs[1::2]=segs_end
                    
        self.segments.extend(tot_segs)

    print "... done loading cells..."
         
    print "...generating cell colours..."
    self.segments_colours = []
    for k in range(len(self.segments)):     #should do this in one step...
        self.segments_colours.append([0,255,0])
    
    print "...converting lists to arrays ..."
    self.segments = np.array(self.segments)
    self.segments_colours = np.array(self.segments_colours)
    print self.segments.shape

    print "... done all ..."


def map_somas(self, cells_select_df):

    self.load_soma_sphere()    #Open sphere primitive

    self.triangle_points=[]
    self.triangle_colours=[]

    for gid, cell_prop in cells_select_df.iterrows():  
        if gid%1000==0: print "... processing cell: ", gid
        
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
        
        soma_start = self.morphologies[model_id]["segs_start"][0]
        soma_end =   self.morphologies[model_id]["segs_end"][0]
        
        soma_start = pos_soma + np.dot(soma_start,RotYZ.T)
        soma_end = pos_soma + np.dot(soma_end,RotYZ.T)
        
       
        size=np.linalg.norm(soma_start-soma_end)          #Size of cell soma;
        
        soma_centre=(soma_end+soma_start)/2.


        #Make triangle surfaces
        for j in range(len(self.triangle_faces)):
            self.triangle_points.append([[self.vertices[int(self.triangle_faces[j][0])-1][0]*size/2.+soma_centre[0],
                                        self.vertices[int(self.triangle_faces[j][0])-1][1]*size/2.+soma_centre[1]-size/2,
                                        self.vertices[int(self.triangle_faces[j][0])-1][2]*size/2.+soma_centre[2]],
                                        
                                         [self.vertices[int(self.triangle_faces[j][1])-1][0]*size/2.+soma_centre[0], 
                                        self.vertices[int(self.triangle_faces[j][1])-1][1]*size/2.+soma_centre[1]-size/2, 
                                        self.vertices[int(self.triangle_faces[j][1])-1][2]*size/2.+soma_centre[2]], 
                                        
                                         [self.vertices[int(self.triangle_faces[j][2])-1][0]*size/2.+soma_centre[0],
                                         self.vertices[int(self.triangle_faces[j][2])-1][1]*size/2.+soma_centre[1]-size/2,
                                         self.vertices[int(self.triangle_faces[j][2])-1][2]*size/2.+soma_centre[2]]                                             
                                         ])

                                
            color_factor = ((self.lowest_vertex_triangle[j])+2.)/4. #Cat: take z direction and normalize
            #self.triangle_colours.append([[int(CMAP[gid%10][0]*color_factor),
            #                                int(CMAP[gid%10][1]*color_factor),
            #                                int(CMAP[gid%10][2]*color_factor)
            #                                ]]*3)
             
            for k in range(3):
                self.triangle_colours.append([int(CMAP[gid%10][0]*color_factor),
                                              int(CMAP[gid%10][1]*color_factor),
                                              int(CMAP[gid%10][2]*color_factor)
                                              ])
    
    
    self.triangle_points = np.array(self.triangle_points) #*1E3
    self.triangle_colours = np.array(self.triangle_colours)


def load_soma_sphere(self):     

    #Load sphere vertices
    f = open('sphere_vertices.csv', 'rt')
    vertices_list = list(csv.reader(f))
    f.close()
    vertices = []
    for row in vertices_list:
        vertices.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3])])
    self.vertices = vertices
    
    #Load sphere faces
    f = open('sphere_faces.csv', 'rt')
    face_list = list(csv.reader(f))
    f.close()
    triangle_faces=[]
    quad_faces=[]
    for row in face_list:
        if (len(row[0].split(" ")[1:]))==4:     #Can use quad spheres or triangle spheres
            quad_faces.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3]), float(row[0].split(" ")[4])])
        else:
            triangle_faces.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3])])

    self.triangle_faces = triangle_faces
    self.quad_faces = quad_faces

    lowest_vertex_triangle = []
    for j in range(len(self.triangle_faces)):
        lowest_vertex_triangle.append(min(self.vertices[int(self.triangle_faces[j][0])-1][1],
                            self.vertices[int(self.triangle_faces[j][1])-1][1],
                            self.vertices[int(self.triangle_faces[j][2])-1][1]))
    self.lowest_vertex_triangle=lowest_vertex_triangle

    #lowest_vertex_quad = []
    #for j in range(len(self.quad_faces)):
        #lowest_vertex_quad.append(min(self.vertices[int(self.quad_faces[j][0])-1][1],
                #self.vertices[int(self.quad_faces[j][1])-1][1],
                #self.vertices[int(self.quad_faces[j][2])-1][1],
                #self.vertices[int(self.quad_faces[j][3])-1][1]))
    #self.lowest_vertex_quad=lowest_vertex_quad
