#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import sip
sip.setapi('QString', 2) #Sets the qt string to native python strings so can be read without weird stuff

from PyQt4 import QtGui 
import sys

from visapi import *

global dendrites_1D, dendrites_1D_colours

dendrites_1D = [] #Initialize to empty lists just to test later; NOT NECESSARY
dendrites_1D_colours = []

def set_morphologies(morphologies):
    global dendrites_1D, dendrites_1D_colours
       
    a = morphologies[morphologies.keys()[0]]['segs_start']
    b = morphologies[morphologies.keys()[0]]['segs_end']
    
    for k in range(len(a)):                                     #Hack way of interleaving two 2D arrays
        dendrites_1D.append(a[k]); dendrites_1D.append(b[k])
    
    dendrites_1D = np.array(dendrites_1D)
    dendrites_1D_colours = dendrites_1D*0. + [0, 255, 0]

    
def set_somata3d():
    pass
    

def initialize(vis):
    app = QtGui.QApplication(sys.argv)
    GUI = GLWindow(vis)
    
    return GUI, app

def run(app):
    sys.exit(app.exec_())
