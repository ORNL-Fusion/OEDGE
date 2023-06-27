Creating a Plasma Grid
======================

OEDGE simulations of the SOL require a field-aligned mesh with which to perform their calculations on. An excellent video demonstrating how the complex curviliner grids used in SOL codes are converted to a simple 2D Cartesian grid can be seen by following `this link <https://drive.google.com/file/d/11c1cVKEBtrhwk9hndMLAhLGfqtpSc8yd/view?usp=sharing>`.  

There are two options for generating a grid for OEDGE:

  1. DG-Carre: This is the most common option. It is the same grid generator used by the 2D SOL fluid code SOLPS-ITER, and so documentation on how to use it is available. 

  2. FUSE: Some time ago a grid generator was built in the IDL programming language to deal with the limited radial extent of DG-Carre grids. It is largely a black box in that it either works or it doesn't, and the grids are specific to OEDGE only. 

Instructions on generating a grid via both methods will be included here. 

DG-Carre Grid
-------------

DG is a GUI used for manipulating and setting up grid making, while Carre is the software that actually makes the grid. Both are very old. To use DG-Carre, we will be setting up your workspace within the SOLPS-ITER directory. This is as standardized a way as there can be in this field, so we will stick to convention.

First, create your solps-iter directory and copy over the needed files:

  .. code-block:: console

    cd /fusion/projects/solps-results
    mkdir [username]
    cd [username]
    mkdir runs post dg-carre
    cd dg-carre
    cp –r /fusion/projects/solps-results/share/dg-carre_share/* .

Next we want to configure your `.bashrc` file to include all the needed aliases. Your `.bashrc` file is in your home directory on iris (`geany ~/.bashrc`). 

  .. code-block:: console

    export OBJECTCODE=”linux.pgf90”
    export TOKAMAK=”D3D”
    alias solpsiter_public='pushd /fusion/projects/codes/solps/SOLPS-ITER/public_code/solps-iter_release_sept2019;source setup.ksh;popd'
    export SOLPSTOP="/fusion/projects/codes/solps/SOLPS-ITER/public_code/solps-iter_release_sept2019"
    export DGHOME="/fusion/projects/solps-results/[username]/dg-carre"
    alias dghome='cd /fusion/projects/solps-results/[username]/dg-carre'
    alias dg='/fusion/projects/solps-results/[username]/dg-carre/DivGeo/dg/linux.ifort64/dg'
    alias e2d= '/fusion/projects/solps-results/[username]/dg-carre/DG/equtrn/e2d'
    alias dg2dg ='/fusion/projects/solps-results/[username]/dg-carre/DG/equtrn/linux.pgf90/dg2dg'
    alias xmatlab='cd /fusion/projects/solps-results/[username]/post/scripts; matlab &'

Your `.bashrc` file is sourced every time you login, but for now since we've modified it we need to source it manually: `source ~/.bashrc`. 

Now we want to to download an equilibrium (i.e., the "gfile"), convert it to a format for DG and then increase the resolution of it.

  .. code-block:: console

    dghome
    cd class
    mkdir d3d_[username]
    cd d3d_[username]
    mkdir 167196_3500
    cd 167196_3500
    idl
    writeg,167196,3500
    exit
    e2d g167196.03500 output.eq
    dg2dg output.eq outputx2.eq
    dg2dg outputx2.eq outputx4.eq
    dg2dg outputx4.eq outputx8.eq

Before we open up DG, we first will want to copy over an .ogr file, which is just a file of the R, Z coordinates of the most limiting surfaces of the DIII-D wall, into the current directory.

  .. code-block::
    cp /fusion/projects/codes/oedge/zamperinis/data/d3d_wall_geometry_june2016.ogr .

  .. note::
    This .ogr file is from June 2016, which is when shot #167196 ran. When doing your own discharge, you will need to use the corresponding wall coordinates as the DIII-D wall frequently changes. Ask around if anyone already has the file made, particularly the SOLPS-ITER people. If one cannot be found, you can extract the coordinates from the gfile. This is perhaps easiest by creating the gfile in the OMFIT EFIT module and copying over the coordinates. In the OMFIT tree, you find the R, Z coordinates at OMFIT["EFIT"]["FILES"]["gEQDSK"]["RLIM"] and "ZLIM". Copy that data over into a text file in the same format as the .ogr file.

Now open up DG. First import the .ogr file: File > Import > Template. Use Ctrl+P to reset the view. Next, convert the template to elements via Commands > Convert > Template to elements. Your screen should look like this:

  .. image:: dg1.png
    :width: 500

Note we have changed the actions on each of the left (L), center (M) and right (R) mouse clicks. It is good practice to make the L and R actions something benign, like Mark and Zoom/Pan, to avoid accidientally performing an action you didn't mean to. As we will learn, DG can be finnicky and so we must work slowly and diligently. 

For our first action, we are changing M to "Reverse normals". Hover over an element and click M to reverse the normal (the direction of the pink line). We need to reverse all the normals, and fortunately there is a shortcut. To perform an action on all elements, hold down shift when clicking the button, e.g., Shift+M. All the lines should be pointing inwards.


