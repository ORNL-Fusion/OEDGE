a=0
b=100
set rmargin 0 
set lmargin 0 
set tmargin 0 
set bmargin 0 

set size square 1,1

set xrange [0:1]
#set yrange[-1.2:1.2]
set y2range[*:*]

plot 'tempsave.txt' u 1:2, -1, 1



pause 1.

a=a+1
if ( a<b ) reread
