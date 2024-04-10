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

# ======================================================================
class fuseTabHub(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        OVERVIEW    = 1
        RESOURCES   = 2        
        CASE_GROUP  = 3
        CASE_OSM    = 4
        CASE_EIRENE = 5
        CASE_DIVIMP = 6
        CASE_SOLPS  = 7

        tab_names = { "Overview":OVERVIEW, "Resources":RESOURCES, "case1":CASE_GROUP, "case2":CASE_EIRENE }

        for key in sorted(tab_names):
            self.label = QLabel(key)

            self.addTab(self.label,key)

#        self.label = QLabel("Shit")
#        self.label.setMinimumSize(200,200)
#        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
#        self.label.setAlignment(Qt.AlignCenter) 

#        self.label2 = QLabel("Piss")

        #self.tab = QTabWidget()

#        self.addTab(self.label ,"Label1")
#        self.addTab(self.label2,"Label2")

