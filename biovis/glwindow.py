#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

#import sip
#sip.setapi('QString', 2) #Sets the qt string to native python strings so can be read without weird stuff

#import numpy as np

#from PyQt4 import QtGui, QtCore     #SIMPLIFY THESE IMPORTS ****************************************************************************

from PyQt4 import QtGui #, QtCore, QtOpenGL   #QtOpenGL installs via "sudo apt-get install python-qt4-gl"
#from OpenGL.GL import *                     #OpenGL installs in ubuntu via "sudo pip install PyOpenGL PyOpenGL_accelerate"
#from OpenGL.GLU import *                    #only used in one location... not clear if necessary

#from simulation import *

from glwidget import GLWidget

#np.set_printoptions(suppress=True)      #Supress scientific notation printing


class GLWindow(QtGui.QWidget):
    
    def __init__(self, fig, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.glWidget = GLWidget(fig, parent=self)
        
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        
        self.setLayout(mainLayout)
        self.setWindowTitle(self.tr("Biovis"))

        self.show()

