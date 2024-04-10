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
class fuseTabHubCase(QTabWidget):

    def __init__(self, num, parent):
#    def __init__(self, num, parent=None):
        super(QTabWidget, self).__init__(parent)

        self.num = num

#        self.parent().parent().parent().parent().statusBar().showMessage("Weird")        
        
        self.initUI()
        
    def initUI(self):

        if self.num == 1:
            
          okButton = QPushButton("OK")
          cancelButton = QPushButton("Cancel")

          hbox = QHBoxLayout()
#          hbox.addStretch(1)
          hbox.addWidget(okButton)
          hbox.addWidget(cancelButton)

          vbox = QVBoxLayout()
#          vbox.addStretch(1)
          vbox.addLayout(hbox)
        
          self.setLayout(vbox)    
        
#        self.setGeometry(300, 300, 300, 150)
#        self.setWindowTitle('Buttons')    
#        self.show()




        
#        dial = QDial()
#        dial.setNotchesVisible(True)
#        spinbox = QSpinBox()
#        layout = QHBoxLayout()
#        layout.addWidget(dial)
#        layout.addWidget(spinbox)
#        self.setLayout(layout)



        
#        self.label = QLabel("Shit")
#        self.label.setMinimumSize(200,200)
#        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
#        self.label.setAlignment(Qt.AlignCenter) 




#        status = self.statusBar()
#        status.showMessage("Hub action", 5000)        
                            
# ======================================================================
class fuseTabHub(QTabWidget):

    def __init__(self, parent):
#    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        MAIN  = 1
        CASE1 = 2
        CASE2 = 3

        tab_names = { " > ":MAIN, "case1":CASE1, "case2":CASE2 }

        for key in sorted(tab_names):

#            self.label = QLabel(key)
            self.tab_hub_case = fuseTabHubCase(1,self)    
            self.addTab(self.tab_hub_case,key)

#        self.label = QLabel("Shit")
#        self.label.setMinimumSize(200,200)
#        self.label.setContextMenuPolicy(Qt.ActionsContextMenu)
#        self.label.setAlignment(Qt.AlignCenter) 

#        self.label2 = QLabel("Piss")

        #self.tab = QTabWidget()

#        self.addTab(self.label ,"Label1")
#        self.addTab(self.label2,"Label2")

