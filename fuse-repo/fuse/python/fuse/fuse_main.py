#!/usr/bin/env python3
#
# from pyth/chap04/calculate.pyw
#
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
#from future_builtins import **

import sys

#sys.path.append("/usr/share/sip/PyQt5")



from math import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *  #(Qt, SIGNAL)
from PyQt5.QtGui  import *  #(QApplication, QDialog, QLineEdit, QTextBrowser,
                            # QVBoxLayout)


                            
# Tab ... things ...
from fuse_experiment import *

from fuse_hub  import *
from fuse_soledge2d  import *
from fuse_jorek  import *
from fuse_transp  import *
from fuse_gpu  import *
from fuse_annova  import *

#import QVTKRenderWindowInteractor

from fuse_cortex import *

# ======================================================================
class fuseTabSimulationMain(QTabWidget):

    def __init__(self, parent):
#    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

# ======================================================================
class fuseTabSimulation(QTabWidget):

    def __init__(self, parent):
#    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        self.main      = fuseTabSimulationMain(self)
        self.hub       = fuseTabHub(self)
        self.soledge2d = fuseTabSoledge2D()
        self.jorek     = fuseTabJorek()
        self.transp    = fuseTabTransp()
        self.gpu       = fuseTabGpu()
        self.annova    = fuseTabAnnova()

        self.addTab(self.main     ," > ")
        self.addTab(self.hub      ,"HUB")
        self.addTab(self.soledge2d,"SOLEDGE2D")
        self.addTab(self.jorek    ,"JOREK")
        self.addTab(self.transp   ,"TRANSP")
        self.addTab(self.gpu      ,"GPU")
        self.addTab(self.annova   ,"ANNOVA")                
    
# ======================================================================
class fuseTabVisualizationMain(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

# ======================================================================
class fuseTabVisualization(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        self.main   = fuseTabVisualizationMain()        
        self.cortex = fuseTabCortex()

        self.addTab(self.main," > ")
        self.addTab(self.cortex,"CORTEX")

# ======================================================================
class fuseTabSetup(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        self.label = QLabel("Shit")
        self.label.setMinimumSize(200,200)
        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.label.setAlignment(Qt.AlignCenter) 

        self.label2 = QLabel("Piss")

        #self.tab = QTabWidget()

        self.addTab(self.label ,"Label1")
        self.addTab(self.label2,"Label2")

# ======================================================================
class fuseTabHelp(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        self.label = QLabel("Shit")
        self.label.setMinimumSize(200,200)
        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.label.setAlignment(Qt.AlignCenter) 

        self.label2 = QLabel("Piss")

        #self.tab = QTabWidget()

        self.addTab(self.label ,"Label1")
        self.addTab(self.label2,"Label2")        

        #status = self.statusBar()    # need to send an event to the main window I think...
        #status.showMessage("Can I help you?", 5000)        

# ======================================================================
class fuseTabMain(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

# ======================================================================
class fuseTabDesign(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

# ======================================================================
class fuseTabUtility(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

# ======================================================================
class fuseTab(QTabWidget):

    def __init__(self, parent):
#    def __init__(self, parent=None):        
        super(QTabWidget, self).__init__(parent)

        self.parent().statusBar().showMessage("Main tabs")        
        
        self.main          = fuseTabMain()
        self.experiment    = fuseTabExperiment()
        self.simulation    = fuseTabSimulation(self)
        self.visualization = fuseTabVisualization()
        self.design        = fuseTabDesign()
        self.utility       = fuseTabUtility()
        self.setup         = fuseTabSetup()        
        self.help          = fuseTabHelp ()

        self.addTab(self.main         ," > ")
        self.addTab(self.experiment   ,"Experiment")
        self.addTab(self.simulation   ,"Simulation")
        self.addTab(self.visualization,"Visualization")
        self.addTab(self.design       ,"Design")
        self.addTab(self.utility      ,"Utility")
        self.addTab(self.setup        ,"Setup")
        self.addTab(self.help         ,"Help")

# ======================================================================
class MainWindow(QMainWindow):
#   --------------------------------------------------------------------
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.setWindowTitle("FUSE")

        self.tab = fuseTab(self)
        #        self.tab = fuseTab()

        self.setCentralWidget(self.tab)

        # Status bar at the bottom of the window:
        #self.sizeLabel = QLabel("Here is the message")
        #self.sizeLabel.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        #status = self.statusBar()
        #status.setSizeGripEnabled(False)
        #status.addPermanentWidget(self.sizeLabel)
#        self.statusBar().showMessage("Loading")#, 5000)        
#   --------------------------------------------------------------------
    #def updateUi(self):
    #    try:
    #        text = unicode(self.lineedit.text())
    #        self.browser.append("{0} = <b>{1}</b>".format(text,
    #                            eval(text)))
    #    except:
    #        self.browser.append("<font color=red>{0} is invalid!</font>"
    #                            .format(text))
#   --------------------------------------------------------------------

# ======================================================================
app = QApplication(sys.argv)
form = MainWindow()
form.show()
app.exec_()

