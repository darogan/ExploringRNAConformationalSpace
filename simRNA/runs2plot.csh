grep '\.' runs.log | awk '{n++; printf("%d %7.3f",m,$1); if (n==3){m++; n=0; print "\n"; if(m==3){m=0}}}' > temp.dat
grep '^0' temp.dat > plot0.dat
grep '^1' temp.dat > plot1.dat
grep '^2' temp.dat > plot2.dat
gnuplot home/runs1.gnu
evince runs.ps
