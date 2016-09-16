set style line  1 lt rgb "black"
set style line  2 lt rgb "red" 
set style line  3 lt rgb "blue" 
set style line  5 lt rgb "green" 
set style line 10 lt rgb "green" 
set style line 31 lt rgb "green"
plot [0:335][0:335] \
'grid.plot' w l ls 31 not, \
'boxA.plot' w l ls 10 lw 2 not, \
'boxP.plot' w l ls 10 lw 2 not, \
'line.plot' u 2:1 w l ls 10 lw 4 not, \
'pairs50.plot'  u 1:2 ls 2  pt 7 ps 1.0 not, \
'pairs200.plot' u 1:2 ls 2  pt 7 ps 0.5 not, \
'pairs50.plot'  u 2:1 ls 2  pt 7 ps 1.0 not, \
'pairs200.plot' u 2:1 ls 2  pt 7 ps 0.5 not, \
'predS.plot'  ls 3 pt 5 ps 1.0 not, \
'preds.plot'  ls 3 pt 5 ps 0.5 not, \
x ls 1 not
set terminal postscript color size 8,7.5
set output 'plot.ps'
replot

