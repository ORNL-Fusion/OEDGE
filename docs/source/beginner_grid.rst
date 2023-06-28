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

    $ dghome
    $ cd class
    $ mkdir d3d_[username]
    $ cd d3d_[username]
    $ mkdir 167196_3500
    $ cd 167196_3500
    $ idl
    $ writeg,167196,3500
    $ exit
    $ e2d g167196.03500 output.eq
    $ dg2dg output.eq outputx2.eq
    $ dg2dg outputx2.eq outputx4.eq
    $ dg2dg outputx4.eq outputx8.eq

Before we open up DG, we first will want to copy over an .ogr file, which is just a file of the R, Z coordinates of the most limiting surfaces of the DIII-D wall, into the current directory.

  .. code-block::

    $ cp /fusion/projects/codes/oedge/zamperinis/data/d3d_wall_geometry_june2016.ogr .

  .. note::
    This .ogr file is from June 2016, which is when shot #167196 ran. When doing your own discharge, you will need to use the corresponding wall coordinates as the DIII-D wall frequently changes. Ask around if anyone already has the file made, particularly the SOLPS-ITER people. If one cannot be found, you can extract the coordinates from the gfile. This is perhaps easiest by creating the gfile in the OMFIT EFIT module and copying over the coordinates. In the OMFIT tree, you find the R, Z coordinates at ``OMFIT["EFIT"]["FILES"]["gEQDSK"]["RLIM"]`` and ``"ZLIM"``. Copy that data over into a text file in the same format as the .ogr file.

Now open up DG with ``dg &``. First import the .ogr file: File > Import > Template. Use Ctrl+P to reset the view. Next, convert the template to elements via Commands > Convert > Template to elements. Your screen should look like something this:

  .. image:: dg1.png
    :width: 500

Note we have changed the actions on each of the left (L), center (M) and right (R) mouse clicks. It is good practice to make the L and R actions something benign, like Mark and Zoom/Pan, to avoid accidentally performing an action you didn't mean to. As we will learn, DG can be finnicky and so we must work slowly and diligently. With that said, save often! File > Save and then name your file something like dg167196.dg.

For our first action, we are changing M to "Reverse normals". Hover over an element and click M to reverse the normal (the direction of the pink line). We need to reverse all the normals, and fortunately there is a shortcut. To perform an action on all elements, hold down shift when clicking the button, e.g., Shift+M. All the lines should now be pointing *away* from the plasma like in the image above.

Next, we import the equlibrium: File > Import > Equlibrium. Select the file outputx8.eq. Now import the topology: File > Import > Topology. The topology options are found in ``$DGHOME/DivGeo/dg/topologies``. If this directory is blank, you may just have a file filter on. At the top bar make sure it ends with a \* and not a filter such as \*.dg (\* is called a wildcard, meaning it will match anything). This discharge uses the SN topology, so select that one. The topology identifies the separatrix. 

Once we've identified the separatrix, we no longer want the red/blue equlibrium on our display as it will quickly get crowded. Let's get a pop-out of the display options as they are handy to have nearby. Click View > Display, and then drag your mouse on the top-most set of dashed lines. This will pop the Display options out to a separate window. Your screen should now look like the following:

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

Before we do the same for Structure, there is one pecularity we must take care of first. DG-Carre requires that the targets be closed polygons. This means we need to create imaginary elements outside of the vessel that make each target a polygon. Whatever shape this polygon is does not matter, just that it needs to be there. To do this, go to Edit > Create > Point. For the outer target, create a point at 1500, -1600. We want to create elements connecting the end points of the outer target to this point (Tip: You can see the end points by Ctrl+U and then clicking Mark in the Structure window for the Outer Target). Set M to Connect Points. Click and drag between the two end points of the outer target. The surface normals must always be facing the inside of the polygon, so switch M to Reverse Normals and click the segments where the normal is facing the wrong way. Your outer target should now look like the following:

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

The Area variable defines the region within which flux surfaces are created. The Cells variable defines the number of flux surfaces for that region. The Delta1 variable defines the spacing between flux surfaces near the X-point. The Delta2 variable defines the spacing between flux surfaces at the other end. You can drag the line to get the desired spacing (see image). Follow the guidelines below when creating flux surfaces:

  - OEDGE can handle finer grid than SOLPS-ITER (for those who have already seen such grids), which generally stick to around 20-30 surfaces in each region. There is no hard and fast rule for OEDGE here, but 20 surfaces within the Core and PFZ and 50 within the SOL is a decent starting point. 
  - Make sure that the number of flux surfaces is the same within the core and the PFZ.
  - Make flux surface spacing small near the X-point and large away from the X-point. The increase in flux surface spacing should be nonlinear.
  - Flux surface spacing should be the same across the X-point in all directions. Users should start with the PFR region and adjust spacing relative to the settings for the PFR region.

  .. image:: dg10.png
    :width: 500

Next is to create grid points, which determine the poloidal resolution of the grid. Go to Edit > Create > Grid Point(s). The Zone variable defines the part of the separatrix along which grid points are created. The Cells variable defines the number of grid points for that zone. The Delta1 variable defines the spacing between grid points near the X-point. The Delta2 variable defines the spacing between grid points at the other end. Make sure that the Law variable is set to “Delta.” You can also drag the plot to get the desired grid spacing. Follow the guidelines below when creating grid points:

  - Keeping the total number of grid points (sum of the outer divertor, inner divertor and SOL) to around 200 is probably a good starting point.
  - Unlike flux surfaces, the number of grid points in each zone do not need to be the same. Make sure that there are a sufficient number of grid points in each divertor zone to yield a high resolution near the target surfaces.
  - Make grid point spacing smaller near the target (Delta2) and larger near the X-point (Delta1) for the divertor zones.
  - Make grid point spacing smaller on both ends for the SOL zone. An easy way to do this is to set Delta1=Delta2. Reduce the Delta value to reduce the grid point spacing near the X-point.
  - The grid point spacing should be the same across the X-point in all directions.

  .. image:: dg11.png
    :width: 500

Now we are just about ready to move onto Carre and to try and generate the actual grid. Before this, go to Commands > Check Variables. This should give a little message in the bottom right of the screen that says, "All variables have valid values". Then click Commands > Rebuild Carre Objects. This has no output, we just trust it does whatever it does correctly. Then click File > Output and press OK. 

Return to your terminal and change to your class directory: ``cd $DGHOME/class``. Run the linking scripting, e.g., ``./lns d3d_[username]/167196_3500/dg167196``. Do not include ``.dg`` in the command. Then navigate to the Carre directory with ``cd $DGHOME/Carre`` and run Carre with ``./carre -``. As an example:

  .. image:: dg12.png
    :width: 500

Carre is run in the following order: Prepare (P), Grid (g), Convert (c), Store (t) and Quit (q). If life was great, we would only need to run this once and we'd have everything we need. Unfortunately, Carre often generates a number of errors and it's largely a black box on how it works (it doesn't help that much of the code is documented in French). The values that work for this tutorial may not necesarilly be the ones that work for the grid you are constructing. Nonetheless, you will encounter this step every time, so it is worth learning how to get around. Some tips:

  - Often decreasing pasmin will solve some errors
  - Pay attention to what zone it has issues with, the error message tells you
  - Often just small changes are all that's needed. So if the error says to change deltr1 and deltrn in zone 2, and you see that they are -4.8467E-4 and -4.213E-2, then try setting them to -0.00001 and -0.001, respectively (i.e., decrease them by a factor of 10). 
  - Sometimes this just comes down to random, dumb luck. I'm sorry.

If you can get through the Grid step, the rest is just hitting Enter at each step until Carre finishes. Next we wish to load the generated grid into DG and inspect it. In DG go to File > Import > Mesh. The grid is stored at $DGHOME/Carre/data/meshes. If you're lucky, you will have the following purple grid loaded into your session.

  .. image:: dg13.png
    :width: 500

Congratulations! You've made your first grid. The file you need for OEDGE will be the one saved as ``$DGHOME/Carre/data/meshes/dg167196.v001.griddivimp``, that is, the ``.griddivimp`` extension is saved in the format needed to be read by OEDGE. You can download the `grid generated by this guide here <https://drive.google.com/file/d/1Egpzb_ZMia2L89t5sh2RjIdl0wmRHvjU/view?usp=sharing>`_. 

Troubleshooting

  - **I loaded my grid into DG, but some of the cells are pink. What does that mean?** Your grid has issues with cells not being orthogonal. Go back into Carre and tweak some of the parameters to generate a new grid. Repeat until the loaded grid has no pink cells. 


FUSE Grid
---------

As an alternative to DG-Carre, OEDGE has support for an additional grid type generated by the scripts within the FUSE repository, collectively called GRID. These are sometimes just refered to as "extended grids" due to the fact GRID grids generally fill in more of the vessel. If one is particularly concerned with plasma far out in the SOL, they may want to consider making the grid with GRID instead of DG-Carre. GRID may also be a bit easier to use since it hides most of the grid making magic behind the scenes, but it is largely a black box as to how it works. In this section we detail how to download the required scripts on iris and how to make the grid.

A memo on GRID can be seen by `clicking this link <https://drive.google.com/file/d/1ElMrd17_plpPB3jyl1tCfgDl3mR7YMSr/view?usp=sharing>`_. We will repeat all the needed instructions for using GRID on this page, but the memo may be useful for anyone needing to run GRID on a machine other than iris. 

Open up your ``.bashrc`` file and add the following lines at the bottom of it:

  .. code-block:: console

    export FUSEHOME=/fusion/projects/codes/oedges/fuse
		alias cdf='cd $FUSEHOME'
		alias cdfo='cd $FUSEHOME/src/osm' 
		alias cdfe='cd $FUSEHOME/src/eirene07' 
		alias cdfi='cd $FUSEHOME/input'
		alias cdfl='cd $FUSEHOME/idl'
		alias cdfs='cd $FUSEHOME/shots'
		alias cdfr='cd $FUSEHOME/results' 
		alias cdfc='cd $FUSEHOME/cases'
		export PATH=$PATH:$FUSEHOME/scripts

    # I have this one commented out...
    #export IDL_PATH=\+$IDL_DIR/lib:+~/fuse/idl:+~/fuse/idl/utils

Next open up your ``.cshrc`` file and add the following at the bottom of it:

  .. code-block:: console

    setenv FUSEHOME "/fusion/projects/codes/oedge/fuse"
    alias cdf  cd $FUSEHOME
    alias cdfo cd $FUSEHOME/src/osm
    alias cdfe cd $FUSEHOME/src/eirene07
    alias cdfi cd $FUSEHOME/input
    alias cdfl cd $FUSEHOME/idl
    alias cdfs cd $FUSEHOME/shots
    alias cdfr cd $FUSEHOME/results
    alias cdfc cd $FUSEHOME/cases
    setenv PATH $PATH":$FUSEHOME/scripts"
    setenv IDL_STARTUP $HOME/idl_startup.pro

Then make sure to source the files after saving it with ``source ~/.bashrc`` and ``source ~/.cshrc``. Now navigate to the fuse directory and create a directory where you will make the grid:

  .. code-block:: console

    $ cdf
    $ cd shots/d3d
    $ mkdir [iris_username]_167196
    $ cd [iris_username]_167196

Next we need to download the gfile into our folder. This can quickly be done with:

  .. code-block:: console

    $ idl
    > writeg,167196,3500
    > exit

There is sometimes a compatability issue with the gfile and GRID depedning on how your download your gfile, but this is not always the case. Sometimes there are not spaces before the mius sigs in the gfile. To address it, we must open the gfile with ``geany g167196.03500 &``. Go to Search > Replace. Make sure "Use regular expressions" is checked. Copy the following regex into the "Search for:" box ``(?<=[0-9])-``. Copy the following into the "Replace with:" box `` -`` (that is a space and then a minus sign). Then click the "In Document" button to add a space before every minus sign.

Now we run the following command to create the needed .equ files:

  .. code-block:: console

    $ fuse -equ d3d [iris_username]_167196 g167196.03500 d3d_167196_3500

You will now see various .equ files within your directory (for those who have made a DG-Carre grid, these are the same .equ files, just likely a bit higher resolution). Now we need to copy over all the needed IDL scripts to run GRID. Run the following command within the ``[iris_username]_167196`` directory:

  .. code-block:: console

    $ cp -r $FUSEHOME/idl .

Next run:

  .. code-block:: console

    $ fuse -make grid-iris [iris_username]_167196

Open up the file ``idl/grid_input.pro`` within your ``[iris_username]_167196`` directory. This is a confusing file, but we only need to change a small section of it. Scroll down to where it says ``'d3d': BEGIN``. This contains the limited amount of input options we control. First, let's specify our wall file. Scroll down to where you see `wall_file` and change the entry to ``wall_file = 'd3d_wall_june2016.dat'``. We will first run the preview option to get a sense of what region our grid may fill in.

  .. code-block:: console

    $ fuse -make grid
    $ fuse -grid-iris -preview d3d [iris_username]_167196 167196_3500.x16.equ test

You should get something like the following output:

  .. image:: grid1.png
    :width: 500

This plot is telling us what the settings in the input file are telling GRID to do. The blue line is setting the radial extent of the generated grid (if it can even be generated that far out without an error). Likewise for the red line in PFZ. Without changing any of the settings, let's see if GRID can succesfully generate a grid. Afterwards we will push the bounds some to get a more extended grid. Run the same command but without ``-preview``:

  .. code-block:: console

    $ fuse -grid-iris d3d [iris_username]_167196 167196_3500.x16.equ test
    
You should have encountered the error ``grid_ANALYSEBOUNDARY``
