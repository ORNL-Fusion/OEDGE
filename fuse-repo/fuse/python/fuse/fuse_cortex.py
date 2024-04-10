#!/usr/bin/env python
#
# from pyth/chap04/calculate.pyw
#


from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
#from future_builtins import **

import sys
from math import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *  #(Qt, SIGNAL)
from PyQt5.QtGui  import *  #(QApplication, QDialog, QLineEdit, QTextBrowser,
                            # QVBoxLayout)

from fuse_cortex_test import *

# ======================================================================
class fuseTabCortexMain(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)
        
# ======================================================================
class fuseTabCortex(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        MAIN  = 1
        CASE1 = 2
        CASE2 = 3
        PLUS  = 4

        tab_names = { " > ":MAIN, "plot1":CASE1, "plot2":CASE2, "+":PLUS }

        tab_type = [1, 2, 3, 4]
        tab_name = [" > ","plot1","plot2","+" ]

        for i in range(4):
#        for key, value in tab_names.items():
#        for key, value in sorted(tab_names.items()):

            self.tab_cortex_plot = fuseTabCortexPlot(tab_type[i])    
            self.addTab(self.tab_cortex_plot,tab_name[i])
#            self.tab_cortex_plot = fuseTabCortexPlot(value)    
#            self.addTab(self.tab_cortex_plot,key)











#        self.label = QLabel("Shit")
#        self.label.setMinimumSize(200,200)
#        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
#        self.label.setAlignment(Qt.AlignCenter) 
#        self.label2 = QLabel("Piss")
        #self.tab = QTabWidget()
#        self.addTab(self.label ,"Label1")
#        self.addTab(self.label2,"Label2")


#        self.main = fuseTabCortexMain()        
#        self.label = QLabel("Shit")
#        self.label.setMinimumSize(200,200)
#        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
#        self.label.setAlignment(Qt.AlignCenter) 
#        self.label2 = QLabel("Piss")
#        #self.tab = QTabWidget()
#        self.addTab(self.main ,">")        
#        self.addTab(self.label ,"Label1")
#        self.addTab(self.label2,"Label2")

