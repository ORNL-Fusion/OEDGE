Beginner's Guide
================

This guide is intended to provide a simple starting point for new users to the OEDGE suite of codes. The codes are designed for use in tokamaks and the generally take full advantage of the simplicity gained by approximating the 3D tokamak (toroidal, radial and poloidal) as a 2D, toroidally-symmetric plasma (poloidal and radial).  By the end of the guide, one should be able to go through the motions of:

  - Generating a computational grid

  - Generate a basic plasma background using OSM-EIRENE

  - Perform a Monte Carlo simulation of impurity transport through a plasma background using DIVIMP

  - Make some basic plots of OSM-EIRENE and DIVIMP output, and learn how to extract the data for futher use

.. note::
  Throughout the literature the naming convention of the code and its different parts is not always consistent. For the purposes of this document, we will use the following:

  - **OEDGE**: The full suite of codes containing all the source code.

  - **OSM-EIRENE**: The 1D fluid solver used to generate a plasma background. Often just abbreviated as OSM.

  - **DIVIMP**: The Monte Carlo trace impurity transport code that follows impurities in a given background.



