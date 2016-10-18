#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import numpy as np

#from PyQt4 import QtGui, QtCore, QtOpenGL   #QtOpenGL installs via "sudo apt-get install python-qt4-gl"
#import OpenGL.GL as GL                      #OpenGL installs in ubuntu via "sudo pip install PyOpenGL PyOpenGL_accelerate"
                                            #Cat: will work on this being explicit rather than wildcard...
#from OpenGL.GLU import gluPerspective       #only used in one location... not clear if necessary

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL, GLU
#from PyQt4.QtOpenGL import QGLShaderProgram

        
class GLWindow(QtGui.QWidget):
    
    def __init__(self, fig, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.glWidget = GLWidget(fig, parent=self)

        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        
        self.setLayout(mainLayout)
        self.setWindowTitle(self.tr("Biovis"))

        self.size = 1

        self.show()

    #def keyPressEvent(self, e):
        
        #if e.key() == QtCore.Qt.Key_Escape:
            #self.close()
            

class GLWidget(QtOpenGL.QGLWidget):
        
    def __init__(self, fig, parent=None):
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.fig = fig

        self.setFocusPolicy(Qt.StrongFocus)

        self.lastPos = QtCore.QPoint()
        self.focus = np.float32([0, 0, 0]) # init camera focus
        self.axes = 'both' # display both mini and focal xyz axes by default
        

        self.load_data_from_fig(fig)
        format = QtOpenGL.QGLFormat()
        
        #format.setVersion(3, 0) # not available in PyQt 4.7.4
        # set to color index mode, unsupported in OpenGL >= 3.1, don't know how to load
        # GL_ARB_compatibility extension, and for now, can't force OpenGL 3.0 mode.
        # Gives "QGLContext::makeCurrent(): Cannot make invalid context current." error:
        #format.setRgba(False)
        
        format.setDoubleBuffer(True) # works fine
        self.setFormat(format)
        #QtOpenGL.QGLFormat.setDefaultFormat(format)
        
        '''
        c = QtGui.qRgb
        cmap = [c(255, 0, 0), c(0, 255, 0), c(0, 0, 255), c(255, 255, 0), c(255, 0, 255)]
        colormap = QtOpenGL.QGLColormap()
        colormap.setEntries(cmap)
        self.setColormap(colormap)
        '''
        
        #self.trigger.connect(self.reloadGL(fig))


    def load_data_from_fig(self, fig):

        self.background = fig.background
        self.segments = np.array(fig.segments)
        self.segments_colours = np.array(fig.segments_colours)
        self.segments3D = np.array(fig.segments3D)
        self.segments3D_colours = np.array(fig.segments3D_colours)
        self.segments3D_joints = np.array(fig.segments3D_joints)
        self.segments3D_joints_colours = np.array(fig.segments3D_joints_colours)
        self.sphere_points = np.array(fig.sphere_points)
        self.sphere_colours = np.array(fig.sphere_colours)
        self.layers = fig.layers
        self.layers_colours = fig.layers_colours
        self.frame = fig.frame
        self.frame_colours = fig.frame_colours

        self.initializeGL()
        self.updateGL()
        self.paintGL()
        
    #def reloadGL(self, fig):
        
        #print "...reloading GL widget..."
        
        #self.segments = fig.segments
        #self.segments_colours = fig.segments_colours
        #self.triangle_points = fig.triangle_points
        #self.triangle_colours = fig.triangle_colours
        #self.layers = fig.layers
        #self.layers_colours = fig.layers_colours
        #self.frame = fig.frame
        #self.frame_colours = fig.frame_colours
        
        #self.initializeGL()
        #self.updateGL()
        #self.paintGL()
        
    def minimumSizeHint(self):
        return QtCore.QSize(50, 50)

    def sizeHint(self):
        return QtCore.QSize(1000, 800)

    def initializeGL(self):
        
        if self.background=='black': GL.glClearColor(0.0, 0.0, 0.0, 1.0) # Black / White toggle switch
        if self.background=='white': GL.glClearColor(1.0, 1.0, 1.0, 1.0)

        GL.glClearDepth(10.0) # same as default
        GL.glEnable(GL.GL_DEPTH_TEST) # display points according to occlusion, not order of plotting
        #GL.glEnable(GL.GL_POINT_SMOOTH) # doesn't seem to work right, proper way to antialiase?
        #GL.glEnable(GL.GL_LINE_SMOOTH) # works better
        #GL.glPointSize(1.5) # truncs to the nearest pixel if antialiasing is off
        #glShadeModel(GL_FLAT)
        #glEnable(GL_CULL_FACE) # only useful for solids
        GL.glTranslate(0, 750, -3000) # init camera distance from origin

        #Modded transparency: blends colours in order they were plotted; otherwise need to use compiled shaders 
        GL.glEnable (GL.GL_BLEND)
        GL.glBlendFunc (GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        #Lighting, not alwasy used
        #GL.glEnable(GL.GL_LIGHTING)
        #GL.glEnable(GL.GL_LIGHT0)
        
        #lightpos = [1000.,0.,0., 1.]
        #GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, lightpos)

        #GL.glDepthMask(GL.GL_TRUE);
        #GL.glDepthFunc(GL.GL_LEQUAL);


        #LightAmbient=[0.2, 0.2, 0.2, 1.0]
        #LightDiffuse=[1.0, 1.0, 1.0, 1.0]
        #LightSpecular= [1.0, 1.0, 1.0, 1.0]
        #LightPosition= [1000.0, 1.0, 1.0]

        #GL.glEnable(GL.GL_LIGHT0)
        #GL.glEnable(GL.GL_LIGHTING)
        #GL.glLightfv(GL.GL_LIGHT0, GL.GL_AMBIENT, LightAmbient)
        #GL.glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE, LightDiffuse)
        #GL.glLightfv(GL.GL_LIGHT0, GL.GL_SPECULAR, LightSpecular)

        #GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, LightPosition)
         
        #GL.glEnable ( GL.GL_COLOR_MATERIAL )
        #GL.glColorMaterial ( GL.GL_FRONT, GL.GL_AMBIENT_AND_DIFFUSE )


        #white = [0.8, 0.8, 0.8, 1.0]
        #cyan = [0., .8, .8, 1.]
        ##GL.glMaterialfv(GL.GL_FRONT, GL.GL_DIFFUSE, cyan)
        #GL.glMaterialfv(GL.GL_FRONT, GL.GL_SPECULAR, white)
        #shininess = [100,100,100,100]
        #GL.glMaterialfv(GL.GL_FRONT, GL.GL_SHININESS, shininess)

        GL.glShadeModel(GL.GL_SMOOTH)

        #self.shaderProgram = QGLShaderProgram()
        #self.shaderProgram.link()


    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        # Don't load identity matrix. Do all transforms in place against current matrix
        # and take advantage of OpenGL's state-machineness.
        # Sure, you might get round-off error over time, but who cares? If you wanna return
        # to a specific focal point on 'f', that's when you really need to first load the
        # identity matrix before doing the transforms
        #GL.glLoadIdentity() # loads identity matrix into top of matrix stack

        # viewing transform for camera: where placed, where pointed, which way is up:
        #GLU.gluLookAt()
        #GL.glScale() # modelling transformation, lets you stretch your objects

        GL.glEnableClientState(GL.GL_COLOR_ARRAY);
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY);
        
        #GL.glEnable(GL_LINE_SMOOTH);
        #GL.glHint(GL_LINE_SMOOTH_HINT,  GL_NICEST);
        


        #Column wire frame
        if len(self.frame)>0:
            GL.glColorPointerub(self.frame_colours) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.frame) # float32
            GL.glDrawArrays(GL.GL_LINES, 0, len(self.frame))



        #Plot regular segments 0-width
        if len(self.segments)>0:
            GL.glColorPointerub(self.segments_colours) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.segments) # float32
            GL.glDrawArrays(GL.GL_LINES, 0, len(self.segments))

            #glColorPointerub(self.dendrite_colors) # unsigned byte, ie uint8
            #glVertexPointerf(self.dendrite_quads) # float32
            #glDrawArrays(GL_QUADS, 0, len(self.dendrite_quads))
           
 
        #Plot 3D segments 
        if len(self.segments3D)>0:
            #GL.glColorPointerub(self.segments3D_colours) # unsigned byte, ie uint8
            #GL.glVertexPointerf(self.segments3D) # float32
            #GL.glDrawArrays(GL.GL_LINES, 0, len(self.segments3D))

            GL.glColorPointerub(self.segments3D_colours) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.segments3D) # float32
            GL.glDrawArrays(GL.GL_QUADS, 0, len(self.segments3D))

            #GL.glColorPointerub(self.segments3D_colours) # unsigned byte, ie uint8
            #GL.glVertexPointerf(self.segments3D) # float32
            #GL.glDrawArrays(GL.GL_QUADS, 0, len(self.segments3D)*4)

        if len(self.segments3D_joints)>0:
            GL.glColorPointerub(self.segments3D_joints_colours) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.segments3D_joints) # float32
            GL.glDrawArrays(GL.GL_QUADS, 0, len(self.segments3D_joints))
            

        #Plots somas 
        if len(self.sphere_points)>0:
            GL.glColorPointerub(self.sphere_colours) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.sphere_points) # float32
            GL.glDrawArrays(GL.GL_TRIANGLES, 0, len(self.sphere_points)*3)  #NOT SURE WHY 3 IS HERE



        #Plot electrode
        if False:
            GL.glColorPointerub(self.colours_electrode) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.electrode) # float32
            GL.glDrawArrays(GL.GL_TRIANGLES, 0, len(self.electrode))
            GL.glColor3ub(255, 255, 255)
            self.renderText (-10*len(self.probe_name),220,0, self.probe_name)


        ##Plots segments or seg3d (3d shapes made out of made of triangles) 
        #if plot_seg3d==1: 
            #GL.glColorPointerub(self.triangle_colours) # unsigned byte, ie uint8
            #GL.glVertexPointerf(self.triangles) # float32
            #GL.glDrawArrays(GL.GL_TRIANGLES, 0, len(self.triangle_points))

        #if movie==1:
        #   GL.glColor3ub(255, 255, 255)
        #   self.renderText (-50,175,0, 'Time: %.1f (ms)' % (float(self.time_index_counter+self.offset)/10.))


        #Plot labels manually
        if False:
            #PLOT Axis tick values
            glColor3ub(255, 255, 255)
            self.renderText (550,0,500, "0.0 um")
            self.renderText (550,-500,500, "-500.0")
            self.renderText (550,-1000,500, "-1000.0")
            self.renderText (550,-1500,500, "-1500.0")
            self.renderText (550,-2000,500, "-2000.0")

            self.renderText (-25,-2060,500, "0.0")
            self.renderText (-530,-2060,500, "-500.0")
            self.renderText (430,-2060,500, "500.0")
                
            if layers_type==0: #plot rat layer indexes
                glColor3ub(255, 255, 255)
                self.renderText (-600,50,500, "Rat")
                self.renderText (-600,-50,500, "L1")
                self.renderText (-600,-175,500, "L2")
                self.renderText (-600,-375,500, "L3")
                self.renderText (-600,-650,500, "L4")
                self.renderText (-600,-925,500, "L5A")
                glColor3ub(0, 127, 255)
                self.renderText (-600,-1200,500, "L5B")
                glColor3ub(255, 255, 255)
                self.renderText (-600,-1500,500, "L6A")
                self.renderText (-600,-1750,500, "L6B")

            if layers_type==1: #plot mouse layer indexes
                self.renderText (-600,50,500, "Mouse")
                glColor3ub(255, 255, 255)
                self.renderText (-600,-50,500, "L1")
                self.renderText (-600,-228,500, "L2/3")
                self.renderText (-600,-414,500, "L4")
                self.renderText (-600,-560,500, "L5")
                self.renderText (-600,-749,500, "L6")


        #text="**********TEST *************"
        #GL.glVertexAttribIPointer(1, 1, GL.GL_UNSIGNED_BYTE, 1, text)
        #GL.glColorPointerub(self.colours_soma) # unsigned byte, ie uint8
        #GL.glDrawArrays(GL.GL_POINTS, 0, len(text))
        
        #print self.points

        if self.axes: # paint xyz axes
            GL.glClear(GL.GL_DEPTH_BUFFER_BIT) # make axes paint on top of data points
            #if self.axes in ['both', 'mini']:
            self.paint_mini_axes()
            #if self.axes in ['both', 'focal']:
            self.paint_focal_axes()

        # might consider using buffer objects for even more speed (less unnecessary vertex
        # data from ram to vram, I think). Apparently, buffer objects don't work with
        # color arrays?



        #Plot layer planes; ************** PLOT LAST for BLENDING ALPHA TRANSPARENCY TO WORK
        if len(self.layers)>0:
            #glColorPointerub(self.layers_colors) # unsigned byte, ie uint8
            #glVertexPointerf(self.layers) # float32
            #glDrawArrays(GL_QUADS, 0, len(self.layers)*4)
            
            
            GL.glColorPointerub(self.layers_colours) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.layers) # float32
            GL.glDrawArrays(GL.GL_TRIANGLES, 0, len(self.layers))


            #glColorPointerub(self.layers_colours) # unsigned byte, ie uint8
            #glVertexPointerf(self.layers) # float32
            #glDrawArrays(GL_TRIANGLES, 0, len(self.layers))



        #GL.glFlush() # forces drawing to begin, only makes difference for client-server?
        self.swapBuffers() # doesn't seem to be necessary, even though I'm in double-buffered
                           # mode with the back buffer for RGB sid encoding, but do it anyway
                           # for completeness

        # print the modelview matrix
        #print(self.MV)

    def save(self, i):
        """Save cluster plot to file"""
        print "SAVING SCREENGRAB"
        fname = "/home/cat/Pictures/1.png"
        
        if fname:
            fname = str(fname) # convert from QString
            image = self.grabFrameBuffer() # defaults to withAlpha=False, makes no difference
            try:
                image.save(fname)
            except Exception as e:
                QtGui.QMessageBox.critical(
                    self.panel, "Error saving file", str(e),
                    QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)

    def resizeGL(self, width, height):
        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        # fov (deg) controls amount of perspective, and as a side effect initial apparent size
        GLU.gluPerspective(45, float(width)/height, 0.1, 100000.) # fov, aspect, nearz & farz
                                                           # clip planes
        #gluPerspective(45, 2.0, 0.1, 100000.) # fov, aspect, nearz & farz
                                                           # clip planes

        GL.glMatrixMode(GL.GL_MODELVIEW)
    

    def paint_mini_axes(self):
        """Paint mini xyz axes in bottom left of widget"""
        w, h = self.width(), self.height()
        vt = self.getTranslation() # this is in eye coordinates
        GL.glViewport(0, 0, w//4, h//4) # mini viewport at bottom left of widget
        self.setTranslation((-1, -.5, -3)) # draw in center of this mini viewport
        self.paint_axes()
        self.setTranslation(vt) # restore translation vector to MV matrix
        GL.glViewport(0, 0, w, h) # restore full viewport

    def paint_focal_axes(self):
        """Paint xyz axes proportional in size to sigma, at focus"""
        GL.glTranslate(*self.focus) # translate to focus
        #self.paint_axes(self.sigma)
        GL.glTranslate(*-self.focus) # translate back

    def update_focal_axes(self):
        """Called every time sigma is changed in main spyke window"""
        #self.update_sigma()
        self.updateGL()

    def paint_axes(self, l=1):
        """Paint axes at origin, with lines of length l"""
        GL.glBegin(GL.GL_LINES)
        GL.glColor3f(1, 0, 0) # red x axis
        GL.glVertex3f(0, 0, 0)
        GL.glVertex3f(l, 0, 0)
        GL.glColor3f(0, 1, 0) # green y axis
        GL.glVertex3f(0, 0, 0)
        GL.glVertex3f(0, l, 0)
        GL.glColor3f(0, 0, 1) # blue z axis
        GL.glVertex3f(0, 0, 0)
        GL.glVertex3f(0, 0, l)
        GL.glEnd()

    def get_MV(self):
        """Return modelview matrix"""
        return GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX) # I think this acts like a copy

    def set_MV(self, MV):
        GL.glLoadMatrixd(MV)

    MV = property(get_MV, set_MV)

    # modelview matrix is column major, so we work on columns instead of rows
    def getViewRight(self):
        """View right vector: 1st col of modelview matrix"""
        return self.MV[:3, 0]

    def getViewUp(self):
        """View up vector: 2nd col of modelview matrix"""
        return self.MV[:3, 1]

    def getViewNormal(self):
        """View normal vector: 3rd col of modelview matrix"""
        return self.MV[:3, 2]

    def getTranslation(self):
        """Translation vector: 4th row of modelview matrix"""
        return self.MV[3, :3]

    def setTranslation(self, vt):
        """Translation vector: 4th row of modelview matrix"""
        MV = self.MV
        MV[3, :3] = vt
        self.MV = MV

    def getDistance(self):
        v = self.getTranslation()
        #return np.sqrt((v**2).sum()) # from data origin
        return np.sqrt(((v-self.focus)**2).sum()) # from focus

    def pan(self, dx, dy):
        """Translate along view right and view up vectors"""
        d = self.getDistance()
        vr = self.getViewRight()
        vr *= dx*d
        GL.glTranslate(vr[0], vr[1], vr[2])
        vu = self.getViewUp()
        vu *= dy*d
        GL.glTranslate(vu[0], vu[1], vu[2])

    def zoom(self, dr):
        """Translate along view normal vector"""
        d = self.getDistance()
        vn = self.getViewNormal()
        vn *= dr*d
        GL.glTranslate(vn[0], vn[1], vn[2])

    def pitch(self, dangle): # aka elevation
        """Rotate around view right vector"""
        vr = self.getViewRight()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vr)
        GL.glTranslate(*-self.focus)

    def yaw(self, dangle): # aka azimuth
        """Rotate around view up vector"""
        vu = self.getViewUp()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vu)
        GL.glTranslate(*-self.focus)

    def roll(self, dangle):
        """Rotate around view normal vector"""
        vn = self.getViewNormal()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vn)
        GL.glTranslate(*-self.focus)

    def panTo(self, p=None):
        """Translate along view right and view up vectors such that data point p is
        centered in the viewport. Not entirely sure why or how this works, figured
        it out using guess and test"""
        if p == None:
            p = self.focus
        MV = self.MV
        vr = self.getViewRight()
        vu = self.getViewUp()
        p = -p
        x = np.dot(p, vr) # dot product
        y = np.dot(p, vu)
        MV[3, :2] = x, y # set first two entries of 4th row to x, y
        self.MV = MV

    def mousePressEvent(self, event):
        self.lastPos = QtCore.QPoint(event.pos())

    def mouseMoveEvent(self, event):
        buttons = event.buttons()
        modifiers = event.modifiers()
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if buttons == QtCore.Qt.LeftButton:
            if modifiers == QtCore.Qt.ControlModifier:
                self.roll(-0.5*dx - 0.5*dy)
            elif modifiers == QtCore.Qt.ShiftModifier:
                self.pan(dx/600., -dy/600.) # qt viewport y axis points down
            else:
                self.yaw(0.5*dx)
                self.pitch(0.5*dy)
        elif buttons == QtCore.Qt.RightButton:
            self.zoom(-dy/500.) # qt viewport y axis points down

        self.updateGL()
        self.lastPos = QtCore.QPoint(event.pos())

    def wheelEvent(self, event):
        self.zoom(event.delta() / 1000.)
        self.updateGL()



