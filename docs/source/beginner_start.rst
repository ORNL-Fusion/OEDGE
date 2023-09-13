Beginner's Guide
================

.. toctree::
  :hidden:

  beginner_workspace
  beginner_grid
  beginner_osm
  beginner_divimp

This guide is intended to provide a simple starting point for new users to the OEDGE suite of codes. The codes are designed for use in tokamaks and the generally take full advantage of the simplicity gained by approximating the 3D tokamak (toroidal, radial and poloidal) as a 2D, toroidally-symmetric plasma (poloidal and radial).  By the end of the guide, one should be able to go through the motions of:

  - Generating a computational grid

  - Generate a basic plasma background using OSM-EIRENE

  - Perform a Monte Carlo simulation of impurity transport through a plasma background using DIVIMP

  - Make some basic plots of OSM-EIRENE and DIVIMP output, and learn how to extract the data for futher use

For now, the tutorial is carried out on the iris cluster at General Atomics and uses DIII-D as an example. You must have access to this machine to go through the guide. Further iterations to the guide could include other machines and devices.

.. note::
  Throughout the literature the naming convention of the code and its different parts is not always consistent. For the purposes of this document, we will use the following:

  - **OEDGE**: The full suite of codes containing all the source code.

  - **OSM-EIRENE**: The 1D fluid solver used to generate a plasma background. Often just abbreviated as OSM.

  - **DIVIMP**: The Monte Carlo trace impurity transport code that follows impurities in a given background.

The general outline of this guide is as follows:

  1. :doc:`beginner_workspace`

  2. :doc:`beginner_grid`

  3. :doc:`beginner_osm`

  4. :doc:`beginner_divimp`

Any questions, comments or concerns about this guide should be directed to Shawn Zamperini (zamperinis@fusion.gat.com). 

