Input Options
=============

OSM and DIVIMP (collectively called OEDGE here) share the same input file as they are two tightly integrated codes that share many of the same options. This inevitably means many options will only apply to OSM and not to DIVIMP, and vice-versa. OEDGE also has a long and rich history of development, thus many input options are considered defunct. 

The input options are separated into "Tags" consisting of a letter and a number, e.g., T13. There are also one set of Tags that start with the number 2 that relate to the SOL22 plasma solver within OSM. Generally, the Tags are supposed to represent a grouping of related input options, but this convention is not always followed. 

`A Tags`_

`A01`_: Title

`A02`_: Case Description

`A03`_: Equilibrium Grid File Name

`200 Tags`_

`201`_: Force Te = Ti

`202`_: Initially Imposed Target Mach Number

A Tags
------

.. _A01:
A01 : Title
  This entry is composed of two string constants. The first is fixed and simply indicates that the first line of the input file is reserved for the title of the run. The second string, after leaving a space, (See the sample input file), contains a title for the case. Typical contents would be the name and series number and perhaps a reference to a note or other document that would describe why the case was run. These are useful in the future when cross-checking the cases that have been run. The title will appear in all print outs and on plots produced from the output of the case.

.. _A02:
A02 : Case Description
  This input is another DIVIMP character string value. The first string is an identifier tag and the second string should contain a description of the important features or options for the case. This description will be included in the header information of the .dat/.html case file and will be included in the description section of the posted case database. 

.. _A03:
A03 : Equilibrium Grid File Name
  This input line specifies the complete path-inclusive name for the equilibrium file that is being used to run the case. This information is then passed to PIN/NIMBUS, if it is being run, for its use in loading the equilibrium grid. At this time, the information is not directly used by DIVIMP itself. DIVIMP relies on the script file to pre-connect the equilibrium grid to the appropriate input unit number.
  e.g. 'Equil File Name' '/u/progs/div4/shots/g31627.v3'

B Tags
------

C Tags
------

D Tags
------

F Tags
------

G Tags
------

H Tags
------

I Tags
------

N Tags
------

P Tags
------

Q Tags
------

S Tags
------

T Tags
------

Z Tags
------

200 Tags
--------

.. _201:
201 : Force Te = Ti
  0: Off - Te and Ti are calculated separately applying the source terms that are appropriate for each species in the independent heat transport equations.

  1: On - Te and Ti are forced to be equal each other at all points - source terms for the two are combined into one heat transport equation.

.. _202 
202 : Initially Imposed Target Mach Number
  This is the value of the flow velocity initially imposed at the target as a multiple of the target sound speed. (A value of 1.0 is usually used initially for the target mach number). If the iterative mach solver option is turned on - then the value of the mach number at the target may move from this initial value as the solver searches for a smooth solution at the point of the super-sonic to sub-sonic transition. 

