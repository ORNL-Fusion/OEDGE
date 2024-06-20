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
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *  #(Qt, SIGNAL)
from PyQt5.QtGui  import *  #(QApplication, QDialog, QLineEdit, QTextBrowser,
                            # QVBoxLayout)


# ======================================================================
class fuseExperimentGrid(QWidget):

    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent)

        dial = QDial()
        dial.setNotchesVisible(True)

        spinbox = QSpinBox()

        layout = QHBoxLayout()
        layout.addWidget(dial)
        layout.addWidget(spinbox)





        self.setLayout(layout)



