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
set log xy
set nolog x
plot 'residus.txt' u :2 w l,"" u :3 w l,"" u :4 w l,"" u :5 w l,"" u :6 w l



pause 1.

a=a+1
if ( a<b ) reread
