import pylab as plt
import numpy as np
import h5py
import math

import matplotlib as mpl
from mpl_toolkits import mplot3d

from matplotlib.collections import LineCollection




def set_axes_equal(ax):
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



def build_collection2d(segs_start,segs_end, color='black'):

    mysegs = []
    
    for seg_start,seg_end in zip(segs_start,segs_end):
        x_start = seg_start[0];   y_start = seg_start[1]   
        x_end = seg_end[0];    y_end = seg_end[1]
    
        mysegs.append([(x_start,y_start),(x_end,y_end)])
    
    line_segments = LineCollection(mysegs,color=color)    

    return line_segments

def build_collections2d(segs_start,segs_end, color='black'):

	'''
	build two collections	
	'''
	xy_segs = []; zy_segs = []

	for seg_start,seg_end in zip(segs_start,segs_end):
		x_start = seg_start[0];   y_start = seg_start[1];   z_start = seg_start[2] 
		x_end = seg_end[0];    y_end = seg_end[1]; z_end = seg_end[2]
		xy_segs.append([(x_start,y_start),(x_end,y_end)])
		zy_segs.append([(z_start,y_start),(z_end,y_end)])


	xy_col = LineCollection(xy_segs,color='k')    
	zy_col = LineCollection(zy_segs,color='k')    

	return [xy_col,zy_col]



def build_collection3d(segs_start,segs_end,color='black'):

	mysegs = []

	for seg_start,seg_end in zip(segs_start,segs_end):
		mysegs.append((seg_start,seg_end))

	segs_col = mplot3d.art3d.Line3DCollection(mysegs,color=color,alpha=0.8)

	return segs_col



def load_morphologies(cm_df,morph_dir):


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




def plot_morph_swc(morph, ax=None):
#	if ax is None:
#		ax = plt.gca()


	ax[0].set_aspect('equal', 'box-forced')
	ax[1].set_aspect('equal', 'box-forced')

	# Make a line drawing of x-y and y-z views
	for p in morph.compartment_list:
		color='black'
		if p['type']==2:
			color='magenta'

		if p['type']==1:
			color='red'

		if p['type']==3:
			color='darkorange'

		
		for c in morph.children_of(p):
			ax[0].plot([p['x'], c['x']], [p['y'], c['y']], color=color)
			ax[1].plot([p['z'], c['z']], [p['y'], c['y']], color=color)

		    
	soma_circle = plt.Circle((morph.soma["x"],  morph.soma["y"]), morph.soma["radius"], color="black", zorder=20)
	ax[0].add_patch(soma_circle)

	soma_circle = plt.Circle((morph.soma["x"],  morph.soma["z"]), morph.soma["radius"], color="black", zorder=20)
	ax[1].add_patch(soma_circle)


	ax[0].set_ylabel('y')
	ax[0].set_xlabel('x')
	ax[1].set_xlabel('z')





def build_collection3d_from_swc(morph,color='black'):

	mysegs = []

	for p in morph.compartment_list:
		for c in morph.children_of(p):
			seg_start = [p['x'],p['y'],p['z']]
			seg_end = [c['x'],c['y'],c['z']]
	
			mysegs.append((seg_start,seg_end))

	segs_col = mplot3d.art3d.Line3DCollection(mysegs,colors=color,alpha=0.8)

	return segs_col





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




def cell_morphs3d(cells_select_df, morphologies,fig,skip=100,color='k'):

	ax = fig.gca(projection='3d')

	col='k'

	RotX = np.array([[1, 0, 0],    # rotate around x axis
				     [0, 0, 1],
				     [0, 1, 0]])

	ncells_select = len(cells_select_df.index)
	rows_display = range(0,ncells_select,skip)
	cells_display_df = cells_select_df.iloc[rows_display]


	for gid,cell_prop in cells_display_df.iterrows():  

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
		pos_soma = np.dot(pos_soma,RotX.T)*1E-3


		segs_start = np.dot(segs_start,RotX.T)*1E-3	
		segs_end = np.dot(segs_end,RotX.T)*1E-3		# for display purposes switch z and y, convert to mm:

		col_segs = build_collection3d(segs_start,segs_end)

		ax.add_collection(col_segs)	# plot segments

		ax.scatter(pos_soma[0], pos_soma[1], pos_soma[2],color='red',alpha=0.4) # plot somata


	set_axes_equal(ax)
	ax.set_aspect('equal')

	ax.set_xlim([-0.8,0.8])
	ax.set_zlim([-1.2,0.2])
	ax.set_ylim([-0.8,0.8])
#	ax.autoscale_view(True,True,True)


	plt.locator_params(nbins=4)

	ax.set_xlabel('lat')
	ax.set_ylabel('lat')
	ax.set_zlabel('depth')




def cell_2d(cell_prop,morphology,ax):

	RotX = np.array([[1, 0, 0],    # rotate around x axis
				     [0, 0, 1],
				     [0, 1, 0]])


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
	pos_soma = np.dot(pos_soma,RotX.T)*1E-3

	segs_start = np.dot(segs_start,RotX.T)	
	segs_end = np.dot(segs_end,RotX.T)		# for display purposes switch z and y, convert to mm:


	xstart = segs_start[:,0]
	ystart = segs_start[:,1]
	zstart = segs_start[:,2]

	xend = segs_end[:,0]
	yend = segs_end[:,1]
	zend = segs_end[:,2]


	LineWidth=0.5
	col='k'

	for iseg in range(len(xstart)):
		if iseg!=0:   
		    ax.plot(plt.r_[xstart[iseg], xend[iseg]], plt.r_[zstart[iseg], zend[iseg]], color=col,lw=LineWidth);plt.hold(True);
		            
		if iseg==0: 
		    ax.scatter(x_soma, y_soma,color='black')





def cell_somata3d(cells_select_df,fig,skip=100,color='yellow'):
	
	ax = fig.gca(projection='3d')

	marker = 'o'
	ncells_select = len(cells_select_df.index)
	rows_display = range(0,ncells_select,skip)
	cells_display_df = cells_select_df.iloc[rows_display]

	x = cells_display_df['x_soma']*1E-3
	y = cells_display_df['y_soma']*1E-3
	z = cells_display_df['z_soma']*1E-3

	ax.scatter(x,z,y,color=color, alpha=0.4,marker=marker)

	set_axes_equal(ax)
	ax.set_aspect('equal')

	ax.set_xlim([-0.8,0.8])
	ax.set_zlim([-1,0.2])
	ax.set_ylim([-0.8,0.8])
#	ax.autoscale_view(True,True,True)

	ax.set_xlabel('x')
	ax.set_ylabel('z')
	ax.set_zlabel('y')


def plot_layer3d(fig):

	import mpl_toolkits.mplot3d as mp3d

	alpha = 0.1

	pia = [(-1, 1, 0),
		   (1, 1, 0),
		   (1, -1, 0),
		   (-1, -1, 0),
		   ]
	l1_bot = [(-1, 1, -0.1),
		   (1, 1, -0.1),
		   (1, -1, -0.1),
		   (-1, -1, -0.1),
		   ]

	l23_bot = [(-1, 1, -0.31),
		   (1, 1, -0.31),
		   (1, -1, -0.31),
		   (-1, -1, -0.31),
		   ]

	l4_bot = [(-1, 1, -0.43),
		   (1, 1, -0.43),
		   (1, -1, -0.43),
		   (-1, -1, -0.43),
		   ]

	l5_bot = [(-1, 1, -0.65),
		   (1, 1, -0.65),
		   (1, -1, -0.65),
		   (-1, -1, -0.65),
		   ]

	l6_bot = [(-1, 1, -0.85),
		   (1, 1, -0.85),
		   (1, -1, -0.85),
		   (-1, -1, -0.85),
		   ]

	planes = [pia,l1_bot,l23_bot,l4_bot,l5_bot,l6_bot]

	ax = fig.gca(projection='3d')
	for plane in planes:
		face = mp3d.art3d.Poly3DCollection([plane], alpha=alpha, linewidth=1)

		# This is the key step to get transparency working
		face.set_facecolor((0, 0, 1, alpha))

		ax.add_collection3d(face)


