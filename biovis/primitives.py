import csv


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


vertices, triangle_faces, lowest_vertex_triangle = load_soma_sphere()

sphere={}
sphere['vertice']=vertices
sphere['triangle_faces']=triangle_faces
sphere['lowest_vertex_triangle']=lowest_vertex_triangle

