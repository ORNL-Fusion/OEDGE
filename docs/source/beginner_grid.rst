Creating a Plasma Grid
======================

OEDGE simulations of the SOL require a field-aligned mesh with which to perform their calculations on. An excellent video demonstrating how the complex curviliner grids used in SOL codes are converted to a simple 2D Cartesian grid can be seen by following `this link <https://drive.google.com/file/d/11c1cVKEBtrhwk9hndMLAhLGfqtpSc8yd/view?usp=sharing>`_.  

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

Next we want to configure your ``.bashrc`` file to include all the needed aliases. Your ``.bashrc`` file is in your home directory on iris (``geany ~/.bashrc``). 

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

Your ``.bashrc`` file is sourced every time you login, but for now since we've modified it we need to source it manually: ``source ~/.bashrc``. 

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
    This .ogr file is from June 2016, which is when shot #167196 ran. When doing your own discharge, you will need to use the corresponding wall coordinates as the DIII-D wall frequently changes. Ask around if anyone already has the file made, particularly the SOLPS-ITER people. If one cannot be found, you can extract the coordinates from the gfile. This is perhaps easiest by creating the gfile in the OMFIT EFIT module and copying over the coordinates. In the OMFIT tree, you find the R, Z coordinates at ``OMFIT["EFIT"]["FILES"]["gEQDSK"]["RLIM"]`` and ``"ZLIM"``. Copy that data over into a text file in the same format as the .ogr file.

Now open up DG. First import the .ogr file: File > Import > Template. Use Ctrl+P to reset the view. Next, convert the template to elements via Commands > Convert > Template to elements. Your screen should look like something this:

  .. image:: dg1.png
    :width: 500

Note we have changed the actions on each of the left (L), center (M) and right (R) mouse clicks. It is good practice to make the L and R actions something benign, like Mark and Zoom/Pan, to avoid accidentally performing an action you didn't mean to. As we will learn, DG can be finnicky and so we must work slowly and diligently. With that said, save often! File > Save and then name your file something like dg167196.dg.

For our first action, we are changing M to "Reverse normals". Hover over an element and click M to reverse the normal (the direction of the pink line). We need to reverse all the normals, and fortunately there is a shortcut. To perform an action on all elements, hold down shift when clicking the button, e.g., Shift+M. All the lines should now be pointing inwards like in the image above.

Next, we import the equlibrium: File > Import > Equlibrium. Select the file outputx8.eq. Now import the topology: File > Import > Topology. The topology options are found in ``$DGHOME/DivGeo/dg/topologies``. If this directory is blank, you may just have a file filter on. At the top bar make sure it ends with a \* and not a filter such as \*.dg (\* is called a wildcard, meaning it will match anything). This discharge uses the SN topology, so select that one. The topology identifies the separatrix. 

Once we've identified the separatrix, we no long want the red/blue equlibrium on our display as it will quickly get crowded. Let's get a pop-out of the display options as they are handy to have nearby. Click View > Display, and then drag your mouse on the top-most set of dashed lines. This will pop the Display options out to a separate window. Your screen should now look like the following:

  .. image:: dg2.png
    :width: 500

Now we must define the target surfaces. These determine how wide your grid will go, at least up until it hits some other limiting surface that prevents the field line from connecting between both targets. Go to Variables > Structure. To set the elements for each target, Mark (we assigned this to L) each segment and then click the Set button in the Structure window. In the below, we have selected the following 3 elements for the inner target (I have turned on the "Points" Display option):

  .. image:: dg3.png
    :width: 500

  .. note::
    Before marking and setting any elements, it is a good habit to always ensure that you have not accidentally marked any other objects. This can create headaches later on down the road. You can unmark all objects, if any are marked, with Ctrl+U. 

For the outer target we set the following elements:

  .. image:: dg4.png
    :width: 500

Before we do the same for Structure, there is one pecularity we must take care of first. DG-Carre requires that the target be closed polygons. This means we need to create imaginary elements outside of the vessel that make each target a polygon. Whatever shape this polygon is does not matter, just that it needs to be there. To do this, go to Edit > Create > Point. For the outer target, create a point at 1500, -1600. We want to create elements connecting the end points of the outer target to this point (Tip: You can see the end points by Ctrl+U and then clicking Mark in the Structure window for the Outer Target). Set M to Connect Points. Click and drag between the two end points of the outer target. The surface normals must always be facing the inside of the polygon, so switch M to Reverse Normals and click the segments where the normal is facing the wrong way. Your outer target should now look like the following:

  .. image:: dg5.png
    :width: 500

Repeat for the inner target by creating another point at 900, -1300 and connecting the points. Pay attention to the normals!

  .. image:: dg6.png
    :width: 500

Now, to set the Structure elements mark ALL elements. Shift+L (assuming L=Mark) allows you to drag a box over everything. Then UNMARK the indicated three elements that separate the targets from the rest of the structure:

  .. image:: dg9.png
    :width: 500

If there were more than just a single element separating the inner and outer targets in the PFZ, you would select just the first elements outside out each defined target, like we are doing for the SOL side of the targets.

Each target needs a SOL and a PFZ edge element designated. Navigate to Variables > Target Specification > #1. Which target corresponds to which target number depends on the configuration:

  - #1: inner target (LSN), outer target (USN), lower inner target (DN)
  - #2: outer target (LSN), inner target (USN), upper inner target (DN)
  - #3: upper outer target (DN)
  - #4: lower outer target (DN)

#167196 is LSN, so Target #1 is the inner target. Mark the furthest out element within the SOL (Ctrl+U first!) and click Set for SOL Edge. Likewise, pick the element furthest into the PFZ and click Set for PFR Edge. Repeat for Target #2 (the outer target). Now is a good reminder to save often if you haven't been!

  .. image:: dg7.png
    :width: 500

We are now ready to create flux surfaces. The first step is to create an innermost flux surface in the core to define the inner extent of the mesh. Set M to Add Surface, and then click somewhere in the core and release to create the inner surface. Then open Edit > Create > Surface(s). 

  .. image:: dg8.png
    :width: 500

The Area variable defines the region within which flux surfaces are created. The Cells variable defines the number of flux surfaces for that region. The Delta1 variable defines the spacing between flux surfaces near the X-point. The Delta2 variable defines the spacing between flux surfaces at the other end. Follow the guidelines below when creating flux surfaces:

  - OEDGE can handle finer grid than SOLPS-ITER (for those who have already seen such grids), which generally stick to around 20-30 surfaces in each region. There is no hard and fast rule for OEDGE here, but 20 surfaces within the Core and 50 within the SOL is a decent starting point. 
  - *SOLPS-ITER only? Make sure that the number of flux surfaces is the same within each region.*
  - Make flux surface spacing small near the X-point and large away from the X-point. The increase in flux surface spacing should be nonlinear.
  - Flux surface spacing should be the same across the X-point in all directions. Users should start with the PFR region and adjust spacing relative to the settings for the PFR region.
  - (Delete) Move the SOL and PFR restraining triangles to adjust the radial extent of the flux
surfaces. However, keep in mind that the flux surfaces must intersect wall elements that
were defined as one of the target surfaces in section 3.2
