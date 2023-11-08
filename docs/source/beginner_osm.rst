Generating an OSM-EIRENE Plasma Background
==========================================

Now we are ready to generate a plasma background using OSM-EIRENE. The process involves building up an input file from nothing, gradually adding various layer of complexity that are required to obtain agreement with experimental data. To begin, we will generate a meaningless background plasma just to demonstrate the workflow.

First `download the grid made by the extended grid generator <https://drive.google.com/file/d/1F3O5wcy5rUo6oAmoXTo5HtM0xLp6pghY/view?usp=sharing>`_ and place it in your ``[iris username]/shots`` directory (e.g., with Filezilla). Using your favorite text editor, such as ``geany``, open up a blank document on iris. We will add a single input option to the input file to indicate we are using a grid from the Fuse grid generator:

  .. code-block:: console
    '+G01    Grid Option      '  3

  .. note::
    **Anatomy of an input option**
    The most basic input options consist of three things: A tag, a description, and a value. In the above, the tag is +G01, the description is "Grid Option" and the value is 3. Every input option has a unique tag and the description is arbitrary and used only to make the input file human-readable. All the input options can be found on this website at :doc:`input`. For instance, documentation for the grid option is found at :ref:`G01`.

Save the input file. The general run command for OEDGE on iris is as follows:

  .. code-block:: console
    $ ./rundiv_master.sh <DIV input file> <OUT input file> <geometry file name> <fluid plasma filename extension - optional> <CFD solution - optional> <DIVIMP solution - optional>"

For our specific instance, we run by replacing the unused files with "none":

  .. code-block:: console
  $ ./rundiv_master.sh d3d-167196-osm-v1 none grid_d3d_167196_3000_v1 none none none

