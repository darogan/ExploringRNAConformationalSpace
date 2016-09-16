set style line  1 lt rgb "green" 
set style line 10 lt rgb "red" 
set style line 19 lt rgb "blue" 
set style line 28 lt rgb "magenta" 
plot [][0:25] \
'plot0.dat' u 2 w l ls  1 lw 3 not, \
'plot0.dat' u 3 w l ls 10 lw 3 not, \
'plot0.dat' u 4 w l ls 19 lw 3 not
set terminal postscript color
set output 'runs.ps'
replot
