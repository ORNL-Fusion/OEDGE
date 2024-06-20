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

# ======================================================================
class fuseTabGpu(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        MAIN   = 1
        CASE1  = 2
        CASE2  = 3

        tab_names = { ">":MAIN, "case1":CASE1, "case2":CASE2 }

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

