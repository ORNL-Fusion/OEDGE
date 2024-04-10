a=0
b=100
set rmargin 0 
set lmargin 0 
set tmargin 0 
set bmargin 0 

set size square 1,1

set xrange [*:*]
set yrange[*:*]
#set y2range[*:*]
set multiplot
#set log xy
#set nolog x
plot 'fluxio3.txt' u :13 w l, "" u :14 w l, "" u :6 w l, "" u :7 w l, "" u :10 w l, "" u :11 w l



pause 1.

a=a+1
if ( a<b ) reread
