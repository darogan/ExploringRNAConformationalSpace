# 1 = cut on interaction (box score), 2 = gap penalty (1.0)
# plot.gnu must exist in home

cat pairs.out  | sort -nr -k 3 | awk '{print $1,$2,"0 8",$3*0.2}' > pairs.dat
home/setsec $argv[1] $argv[2] > setsec.log
sort -n -k3 line.dat | sed 's/2$/2#/' | tr "#" "\n" > line.plot
cp pairs.dat pairs.plot
head -10 pairs.plot > pairs50.plot
head -100 pairs.plot > pairs200.plot
grep ' [67]$' pred.dat > predS.plot
grep -v ' [67]$' pred.dat > preds.plot
# eval "sed 's/NNN/$len/g' ../plot.gnu > plot.gnu"
gnuplot plot.gnu
evince plot.ps &
awk '{if($4>$5){print $3,($4+1)/($5+1)}else{print $3,($5+1)/($4+1)}}' pack.dat > pack.plot
