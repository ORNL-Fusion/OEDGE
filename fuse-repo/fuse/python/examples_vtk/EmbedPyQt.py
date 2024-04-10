#!/usr/bin/env python
 
import sys
import vtk
#from PyQt5 import QtCore, QtGui

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *  #(Qt, SIGNAL)
from PyQt5.QtGui  import *  #(QApplication, QDialog, QLineEdit, QTextBrowser,
                            # QVBoxLayout)

#from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import QVTKRenderWindowInteractor
from QVTKRenderWindowInteractor import *

class MainWindow(QMainWindow):
 
    def __init__(self, parent = None):
        QMainWindow.__init__(self, parent)
 
        self.frame = QFrame()
 
        self.vl = QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)
 
        self.ren = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
 
        # Create source
        source = vtk.vtkSphereSource()
        source.SetCenter(0, 0, 0)
        source.SetRadius(5.0)
 
        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
 
        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
 
        self.ren.AddActor(actor)
 
        self.ren.ResetCamera()
 
        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)
 
        self.show()
        self.iren.Initialize()
 
 
if __name__ == "__main__":
 
    app = QApplication(sys.argv)
 
    window = MainWindow()
 
    sys.exit(app.exec_())
