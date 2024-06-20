#!/usr/bin/env python
#
# from pyth/chap04/calculate.pyw
#
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from future_builtins import *

import sys
from math import *
from PyQt4.QtCore import *  #(Qt, SIGNAL)
from PyQt4.QtGui  import *  #(QApplication, QDialog, QLineEdit, QTextBrowser,
                            # QVBoxLayout)

# Tab ... things ...
from fuse_tab_hub  import *
from fuse_tab_data import *

# ======================================================================
class fuseTabCodes(QTabWidget):

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

        self.hub   = fuseTabHub  ()
        self.data  = fuseTabData ()
        self.codes = fuseTabCodes()
        self.setup = fuseTabSetup()
        self.help  = fuseTabHelp ()

        self.addTab(self.hub  ,"Hub"  )
        self.addTab(self.data ,"Data" )
        self.addTab(self.codes,"Codes")
        self.addTab(self.setup,"Setup")
        self.addTab(self.help ,"Help" )

# ======================================================================
class MainWindow(QMainWindow):
#   --------------------------------------------------------------------
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.setWindowTitle("FUSE")

        self.tab = fuseTabMain()

        self.setCentralWidget(self.tab)

        # Status bar at the bottom of the window:
        #self.sizeLabel = QLabel("Here is the message")
        #self.sizeLabel.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        status = self.statusBar()
        #status.setSizeGripEnabled(False)
        #status.addPermanentWidget(self.sizeLabel)
        status.showMessage("Loading", 5000)        
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

