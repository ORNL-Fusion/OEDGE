# OEDGE Documentation

https://oedge.readthedocs.io/en/latest/index.html

Here we will describe the OEDGE documentation hosted via ReadTheDocs and how one can add/modify it. The documentation is compiled via ReadTheDocs, which uses the sphinx python package to generate webpages. The relevant files are:

- `requirements.txt` - File containing the sphinx package as well as the readthedocs theme needed for sphinx. A specific sphinx version can be requested for stability purposes. The current setting is to use the newest version to avoid falling out of compliance. 
- `source/index.rst` - This acts as the homepage for the documentation. It should include links to all "header" type of pages. E.g., there should be a link to the start of the Beginner's Guide, but a link to every ssection of the guide is not necessary.
- `source/conf.py` - Contains configuration settings for the documentation. Generally nothing should need changes here unless the overall theme of the guide is to be changed.
- `source/*.rst` - The .rst files are the individual pages for the documentation.

## Layout and adding new pages

Pages are linked by designating a page as a starting point, or "header", for a particular section of the website. Consider the layout of the Beginner's Guide. We designated the start page for the guide as `beginner_start.rst` and include it in the table of contents (`toctree`) within `index.rst`:

```
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   beginner_start
   common_workflows
   python_plotting
   input
   publications
```

Then within `beginner_start.rst` we include links to rest of the Beginner's Guide:

```
Beginner's Guide
================

.. toctree::
  :hidden:

  beginner_workspace
  beginner_grid
  beginner_osm
  beginner_divimp
```

In this way a tree of navigable pages is created. When creating new pages this structure should be followed.

## Modifying pages

Changes to the documentation can either be made directly in the GitHub text editor (recommended) or treated lke any other file within the repository via the usual git version control tools. Changes within the GitHub text editor are possible only if you have write access to the branch you are on. This most likely means you will need to make changes on a separate branch and then create a pull request to merge the changes into the master branch. This can all be done on GitHub.

To modify or add to a page, find the corresponding `.rst` file and make the changes accordingly. Most of the relevant syntax for the page can be found here:

https://sphinx-rtd-theme.readthedocs.io/en/stable/demo/demo.html

If making changes within the GitHub, make sure you are on a branch in which you have write access and then press the green Commit Changes... button at the top of the editor. Then click the Pull Requests tab and request to merge your branch into master with the updates. Make sure to add a brief description of the changes. 

## Recompiling the documentation

To recompile the documentation you must have a ReadTheDocs account and be added as a contributor to the OEDGE project. Contact one of the OEDGE developers if you would like to be added, otherwise you can simply notify one of the OEDGE developers and request they recompile the documentation. 

Within the OEDGE project directory on ReadTheDocs, in the Overview section, is a Build Project button. Click this and the proiject will automatically recompile. If there are no errors, then the process is finished and any changes are live.
