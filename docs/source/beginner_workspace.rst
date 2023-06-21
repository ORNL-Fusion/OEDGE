Setting up Your Workspaces
==========================

As mentioned before, this guide sets a user up on the iris cluster at General Atomics. This is because different clusters require different permissions, and it is therefore impossible to make a one-guide-fits-all. Nonetheless, the procedure here is general and it may be possible to read through and adapt for different environments. 

Local Machine Setup
-------------------

First, we will clone the OEDGE repository to your own local machine (e.g., your laptop). The code will not be run on your machine, but the python scripts that are used to interface with the output are within the repository.

Setup on iris
-------------

To begin, you must first be added to the 'oedge' user group on iris. Contact the server administrators (fus@fusion.gat.com) and request to a) be added to the 'oedge' group and b) request access to the ADAS database. You will be sent a User Agreement form for the ADAS data. Once you are added to the 'oedge' group, navigate to ``/fusion/projects/codes/oedge`` and create the following directories:

  1. ``mkdir [iris username]``

  2. ``cd [iris username]``

  3. ``mkdir cases shots results data``

This directory is where all your OEDGE work will live. 

Next you will want to copy over the scripts required to run the codes within OEDGE. Within your OEDGE directory (``/fusion/projects/codes/oedge/[iris username]``) copy:

  1. ``cp /fusion/projects/codes/oedge/zamperinis/rundiv_master.iris .``

  2. ``cp /fusion/projects/codes/oedge/zamperinis/rundiv.sh .``

The first script is the main shell script that runs either OSM-EIRENE or DIVIMP. Since iris uses a scheduler to submit jobs, like most clusters, a slurm script is needed to submit jobs to the queue. That's what the second script is for. 

.. note::

  OSM-EIRENE and DIVIMP share a lot of the same code. In fact, the two codes share the same input file, that's how interconnected they are! This is why you may often hear a background was generated in DIVIMP, but in reality it was made using OSM-EIRENE. 

Now we need to tell the scripts to run in your directory. Open up both files in your text editor of choice (I prefer ``geany``) and change the following lines:

  - In ``rundiv_master.iris``:
    
    - ``PROGDIR=master``
    - ``DATAROOT=/fusion/projects/codes/oedge/[iris username]``
    - ``RUNROOT=/fusion/projects/codes/oedge/[iris username]``
    - ``RESULTSROOT=/fusion/projects/codes/oedge/[iris username]``

  - In ``rundiv.sh``:
   
    - ``echo "srun ./rundiv_master.iris $1 $4 $2 $5 none $3" >> $1.sh``   (we're adding the ``_master``)
