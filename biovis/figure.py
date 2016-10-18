#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import sys
from glwindow import GLWindow
from PyQt4 import QtCore, QtGui

import netgraph as ng       #Sergey's dataframe loading modules
import artist as art    #functions to load data into opengl-ready arrays

class Figure(object):
    '''
	Main API functions
    '''

    def __init__(self):
        print "... initializing canvas ..."
        #self.app = QtGui.QApplication(sys.argv)

        self.app = QtCore.QCoreApplication.instance()
        if self.app is None:
            self.app = QtGui.QApplication([])
    
        self.app.processEvents()
        
        #Set default parameters if user does not set these values; #opengl needs defaults for plotting routines; 
        self.set_defaults()
    
    def set_defaults(self):
        
        ''' Functions to initialize arrays; may not be necessary. 
        '''
        
        self.background = 'black'
        self.segments = []
        self.segments_colours = []
        self.segments3D = []
        self.segments3D_colours = []
        self.segments3D_joints = [] 
        self.segments3D_joints_colours = []        
        self.sphere_points = []
        self.sphere_colours = []
        self.layers = []
        self.layers_colours = []
        self.frame = []
        self.frame_colours = []


    def set_bgcolor(self, color):
        print "...setting gbcolor: ", color
        self.background = color


    def set_frame(self, frame, frame_colours):
        print "...setting frame..."

        if frame == []: self.frame = []; self.frame_colours = []                #Reset somas
        else: self.frame, self.frame_colours = art.draw_frame(frame, frame_colours)
                        

    def set_layers(self, layer_depths, layer_colours, layer_alpha):
        print "...setting layers..."

        if layer_depths == []: self.layers = []; self.layers_colours = []       #Reset morphs
        else: self.layers, self.layers_colours = art.draw_layers(layer_depths, layer_colours, layer_alpha)
       

    def plot_somas(self, cells_select_df, morphologies, cmap,color_label):

        sphere_points, sphere_colours = art.draw_somas(cells_select_df, morphologies, cmap, color_label)
        self.sphere_points.extend(sphere_points)
        self.sphere_colours.extend(sphere_colours)

 
    def plot_morph(self, cells_select_df, morphologies, cmap, color_label):
        
        segments, segments_colours = art.draw_morphologies(cells_select_df, morphologies,cmap,color_label)
        self.segments.extend(segments)
        self.segments_colours.extend(segments_colours)


    def plot_morph3D(self, cells_select_df, morphologies, cmap, color_label, n_faces):
        
        segments3D, segments3D_colours, segments3D_joints, segments3D_joints_colours = art.draw_morphologies3D(cells_select_df, morphologies,cmap,color_label, n_faces)
        self.segments3D.extend(segments3D)
        self.segments3D_colours.extend(segments3D_colours)
        self.segments3D_joints.extend(segments3D_joints)
        self.segments3D_joints_colours.extend(segments3D_joints_colours)


    def plot_slice(self, cells_select_df, morphologies, cmap, color_label, xplane_range):
        
        segments, segments_colours = art.draw_slice(cells_select_df, morphologies,cmap,color_label,xplane_range)
        self.segments.extend(segments)
        self.segments_colours.extend(segments_colours)
    
    def plot_synapses(self, cid, synapses, cells_select_df, morphologies, cmap, color_label, n_faces):
    
        #Plot synapses along surface of cylinders; 
        sphere_points, sphere_colours = art.draw_synapses(self.segments3D, cid, synapses, cells_select_df, morphologies, cmap, color_label, n_faces)
        self.sphere_points.extend(sphere_points)
        self.sphere_colours.extend(sphere_colours)
        

    def clear(self):
        self.set_defaults()


    def set_title(self, title):
        print title
        pass


    def resize_screen(self):
        
        sizes = [1000, 10000]
        
        self.GUI.resize(sizes[self.GUI.size%2], sizes[self.GUI.size%2])
        self.GUI.size+=1

        self.GUI.glWidget.updateGL()

    def screen_grab(self):
        
        self.GUI.glWidget.save(0)

    def show(self):   #Show
        
        print "... showing ..."

        #self.GUI = GLWindow(self)   #Pass entire figure object in order to access its attributes inside opengl
        #self.app.exec_()
        
        self.GUI = GLWindow(self)
        #GUI.show()
        
        #try:
        #from IPython.lib.guisupport import start_event_loop_qt4
        #start_event_loop_qt4(self.app)
        
        #except ImportError:
        self.app.exec_()
        
        #sys.exit()


        #app = QApplication(sys.argv)
        #widget = GLWindow(self)
        #widget.show()
        #sys.exit(self.app.exec_())


                
    def update(self):   #Restarts widget

        #self.GUI.glWidget.repaint()
        #self.GUI.show()
        self.GUI.glWidget.updateGL()
        #self.GUI.glWidget.updateGL()
        #self.GUI.glWidget.repaint()
        #self.GUI.glWidget.show()

        #self.app.processEvents()

        #print "...updating ... NOT WORKING "

        #self.GUI = GLWindow(self)
        #self.app.exec_()
        
