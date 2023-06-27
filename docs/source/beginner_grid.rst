Creating a Plasma Grid
======================

OEDGE simulations of the SOL require a field-aligned mesh with which to perform their calculations on. An excellent video demonstrating how the complex curviliner grids used in SOL codes are converted to a simple 2D Cartesian grid can be seen by following `this link <https://drive.google.com/file/d/11c1cVKEBtrhwk9hndMLAhLGfqtpSc8yd/view?usp=sharing>`.  

There are two options for generating a grid for OEDGE:

  1. DG-Carre: This is the most common option. It is the same grid generator used by the 2D SOL fluid code SOLPS-ITER, and so documentation on how to use it is available. 

  2. FUSE: Some time ago a grid generator was built in the IDL programming language to deal with the limited radial extent of DG-Carre grids. It is largely a black box in that it either works or it doesn't, and the grids are specific to OEDGE only. 

Instructions on generating a grid via both methods will be included here. 

DG-Carre Grid
-------------

