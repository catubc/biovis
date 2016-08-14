#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import sip
sip.setapi('QString', 2) #Sets the qt string to native python strings so can be read without weird stuff

import sys
from PyQt4 import QtGui, QtCore     #MINIMIZE THESE IMPORTS ****************************************************************************
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import numpy as np
from simulation import *

np.set_printoptions(suppress=True)      #Supress scientific notation printing

#****************************************************************************

class Load(QtGui.QWidget):
    def __init__(self, parent):
        super(Load, self).__init__(parent)
        layout = QtGui.QGridLayout()
        self.parent=parent
        
    def load_sim(self, main_widget):
        print("....loading sim ...")
        self.parent.root_dir = 'simulations/' 
        self.simulation = Simulation(QtGui.QFileDialog.getExistingDirectory(main_widget, 'Load Experiment', self.parent.root_dir), self.parent.root_dir)
        self.parent.setWindowTitle(self.simulation.name)


class View_sim(QtGui.QWidget):        
    def __init__(self, parent):
        super(View_sim, self).__init__(parent)
        self.parent=parent
        layout = QtGui.QGridLayout()

        row_index = 0
        #*****************************************************************************
        #************************** SELECT SIMULATION  *******************************
        #*****************************************************************************
        self.select_lbl = QLabel('Select Simulation----->', self)
        self.select_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        layout.addWidget(self.select_lbl, row_index,0)
        self.comboBox_select_simulation = QtGui.QComboBox(self)
        file_names = sorted(glob.glob(self.parent.root_dir+"simulations/*"))
        
        for file_name in file_names:
            self.comboBox_select_simulation.addItem(file_name.replace(self.parent.root_dir+"simulations/",''))
        layout.addWidget(self.comboBox_select_simulation, row_index,1)
        self.comboBox_select_simulation.activated[str].connect(self.select_simulation); self.selected_simulation = file_names[0].replace(self.parent.root_dir+"simulations/",'')
        

        row_index+=1
        

        #*****************************************************************************
        #***************************** PARAMETERS ************************************
        #*****************************************************************************
        
        self.rarifier = QLineEdit("1")
        self.rarifier_lbl = QLabel('% Cells to View', self)
        layout.addWidget(self.rarifier_lbl, row_index, 0)
        layout.addWidget(self.rarifier, row_index,1)
        
        
        row_index+=1


        #*****************************************************************************
        #***************************** QUERY *****************************************
        #*****************************************************************************
        
        self.query_string = QLineEdit("ei=='excitatory'")
        self.query_button = QPushButton('Load Query')
        self.query_button.setMaximumWidth(200)
        self.query_button.clicked.connect(self.ld_query)
        layout.addWidget(self.query_button, row_index, 0)
        layout.addWidget(self.query_string, row_index,1)
        
        
        row_index+=1

        
        #*****************************************************************************
        #***************************** VIEW *****************************************
        #*****************************************************************************
        
        self.view_button = QPushButton('View Network')
        self.view_button.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.view_button.setStyleSheet('color: blue')
        
        self.view_button.clicked.connect(self.vw_network)
        layout.addWidget(self.view_button, row_index, 0)
        
        
        self.setLayout(layout)



        
        #Load default simulation
        self.parent.root_dir = 'simulations/' 
        self.parent.simulation = Simulation(self.selected_simulation, self.parent.root_dir)
        self.ld_query()
        
        self.parent.setWindowTitle(self.parent.simulation.name)
        
        

    def select_simulation(self, text):
        self.parent.root_dir = 'simulations/' 
        self.simulation = Simulation(text, self.parent.root_dir)
        self.parent.setWindowTitle(self.simulation.name)
        
        
    def ld_query(self):
        self.parent.simulation.load_query(self)

    
    def vw_network(self):
        self.parent.simulation.load_query(self)
        self.parent.simulation.show_query()
        
        
        
class FileDialog(QtGui.QFileDialog):
    """ Hack; overwrited intrinsic FileDialog functions to load multiple directories."""
    def __init__(self, *args, **kwargs):
        super(FileDialog, self).__init__(*args, **kwargs)
        self.setOption(QtGui.QFileDialog.DontUseNativeDialog, True)
        self.setFileMode(QtGui.QFileDialog.ExistingFiles)
        #self.setFileMode(QtGui.QFileDialog.DirectoryOnly)

        self.out_files=[]

    def accept(self):
        self.out_files = self.selectedFiles()   #can't index into this directly, ok to use "out_files"
        self.deleteLater()
        #super(FileDialog, self).accept()
        
        
            
class Window(QtGui.QMainWindow):

    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        #self.setWindowTitle("OpenNeuron")
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))


        #Set widget to show up with viewbox
        toolMenu = QtGui.QMenuBar()
        toolMenu.setNativeMenuBar(False) # <--Sets the menu within the widget; otherwise shows up as global (i.e. at top desktop screen)
        self.setMenuBar(toolMenu)

        #***** TEXT PARAMETERS FIELDS ******
        #Mouse July 11 as default experiment
        self.root_dir = '' 

        #Menu Item Lists
        self.make_menu()

        #LOAD CENTRAL WIDGET
        self.central_widget = QtGui.QStackedWidget()
        self.setCentralWidget(self.central_widget)
        
        #SET DEFAULT WIDGET TO PROCESS
        self.view_widget = View_sim(self)
        self.central_widget.addWidget(self.view_widget)
        
        
        self.show()
   



    def make_menu(self):
        
        #DEFINE MENUS
        loadSimulation = QtGui.QAction("&Load Simulation", self)
        loadSimulation.setStatusTip('Load Simulation')
        loadSimulation.triggered.connect(self.ld_simulation)

        exitApplication = QtGui.QAction("&Exit Application", self)
        exitApplication.setStatusTip('Exit')
        exitApplication.triggered.connect(self.close_application)
                
        viewSimulation = QtGui.QAction("&View Simulation", self)
        viewSimulation.setStatusTip('View Simulation')
        viewSimulation.triggered.connect(self.view_simulation)

        
        #MAKE MENUS
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(loadSimulation)
        fileMenu.addAction(exitApplication)
        fileMenu = mainMenu.addMenu('View')


    #************* LOAD FILE MENUS *****************
    def ld_simulation(self):
        self.load_widget = Load(self)
        self.load_widget.load_sim(self)   #Pass main widget to subwidgets as it contains needed parameters.

        #Send to view widget
        self.view_widget = View_sim(self)
        self.central_widget.addWidget(self.view_widget)
        self.central_widget.setCurrentWidget(self.view_widget)


    def view_simulation(self):
        self.load_widget = Load(self)
        self.load_widget.load_sim(self)   #Pass main widget to subwidgets as it contains needed parameters.
        
        self.view_widget = View_sim(self)
        #self.central_widget.addWidget(self.view_widget)
        self.central_widget.setCurrentWidget(self.view_widget)

        
    def close_application(self):
        print("Do svidaniya")
        
        sys.exit()



def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())


run()
