#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import sys
from PyQt4.QtGui import QApplication
from glwindow import GLWindow

import netgraph as ng       #Sergey's dataframe loading modules

class Figure(object):
    '''
	Main API functions
    '''

    def __init__(self):
        print "... initializing canvas ..."
        self.app = QApplication(sys.argv)

        #Set default parameters if user does not set these values; #opengl needs defaults for plotting routines; 
        self.set_defaults()
    
    def set_defaults(self):
        
        self.background = 'black'
        self.segments = []
        self.segments_colours = []
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
        else: self.frame, self.frame_colours = ng.load_frame(frame, frame_colours)
                        

    def set_layers(self, layer_depths, layer_colours, layer_alpha):
        print "...setting layers..."

        if layer_depths == []: self.layers = []; self.layers_colours = []       #Reset morphs
        else: self.layers, self.layers_colours = ng.load_layers(layer_depths, layer_colours, layer_alpha)


    def plot_somas(self, cells_select_df, morphologies):

        self.sphere_points, self.sphere_colours = ng.plot_somas(cells_select_df, morphologies)
 
 
    def plot_morph(self, cells_select_df, morphologies):
        
        self.segments, self.segments_colours = ng.plot_morphologies(cells_select_df, morphologies)
        

    def plot_presynaptic_somas(self, cells_select_df, morphologies):
        
        self.sphere_points, self.sphere_colours = ng.plot_somas(cells_select_df, morphologies)
        
        

    def clear(self):
        self.set_defaults()


    def set_title(self, title):
        print title
        pass


    def show(self):   #Show
        
        self.GUI = GLWindow(self)   #Pass entire figure object in order to access its attributes inside opengl

        print "... showing ..."
        self.app.exec_()

                
    def restart(self):   #Restarts widget
        print "...restarting ..."
        self.GUI = GLWindow(self)
        self.app.exec_()


