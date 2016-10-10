import numpy as np
import math

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



def compute_segs(morphology,cell_prop):

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
    morph_segs_start = morphology["segs_start"]
    morph_segs_end =   morphology["segs_end"] 
        
    segs_start = pos_soma + np.dot(morph_segs_start,RotYZ.T)  # use the tranposed matrix since we multiply on the right
    segs_end = pos_soma + np.dot(morph_segs_end,RotYZ.T)

    return segs_start,segs_end
