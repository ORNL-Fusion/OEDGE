$
$ 200 and 210 series plots are LOS plots from a given 
$ observation position for density and PLRPS respectively.
$ They require a second line of plot data containing:
$ Robs , Zobs , Thetamin(degs), Thetamax(degs), >
$ Drad - the radial integration size (dr in meters) 
$
$ Numthe - number of angular bins for the plot
$
$ Numav  - the number of arcs to average over for a bin value
$
$ Nsmooth- the number of bins to average over for smoothing
$          0 = off
$
$ IZmin,IZmax - range of plots to produce from -2 to nizs+1
$
$ All of the above with an opening character string appear
$ on the same line. The second data line should be commented 
$ out when the plot is turned off. Any value not equal to 0 
$ on the first line will turn the plot on.
$
'200 Density LOS plot                       '     0
'000 Data'  4.0  0.0  75.0  115.0  0.001  100 1 0  0  6
'210 PLRP LOS plot                          '     0
'000 Data'  4.0  0.0  75.0  115.0  0.001  100 1 0  0  6



