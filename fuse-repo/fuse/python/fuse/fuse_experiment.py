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

from fuse_grid import *

# ======================================================================
class fuseTabExperiment(QTabWidget):

    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)

        GRID             = 1
        DIVERTOR_PROBES  = 2        
        UPSTREAM_THOMSON = 3
        DIVERTOR_THOMSON = 4
        CAMERA_2D        = 5

        tab_names = { "Grid":GRID, "LP":DIVERTOR_PROBES, "Thomson":UPSTREAM_THOMSON, "Cameras":CAMERA_2D }

        for key in sorted(tab_names):
            self.grid = fuseGrid()

            self.addTab(self.grid,key)




