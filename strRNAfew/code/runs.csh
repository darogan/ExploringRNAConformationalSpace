# 1 = number of models

#@ n = `cat pairs.pred.dat | wc -l`

rm -r save*
rm pairs.log

cat code/stem.model | sed 's/MOVE/   10/' > stem.model
@ len = `cat true.pdb | wc -l`
@ n = $len / 3 + 10
#@ n = $len / 4 + 5
#@ n = ( 10 * $len ) / 35
echo using $n pairs
echo using $n pairs > runs.log
echo >> runs.log
@ pairs = $n + 1
while ( $pairs > 0 )
	@ pairs--
	echo pairs = $pairs
	echo >> runs.log
	echo pairs = $pairs >> runs.log
	tcsh code/pairs.csh $argv[1] $pairs >> runs.log
	set grem = `cat pairs.grem.rms | awk '{n++; s+=$3; print s/n}' | tail -1`
	set pred = `cat pairs.pred.rms | awk '{n++; s+=$3; print s/n}' | tail -1`
	set true = `cat pairs.true.rms | awk '{n++; s+=$3; print s/n}' | tail -1`
	echo $pairs $grem $pred $true >> pairs.log
	mkdir save$pairs
	cp models.*.pdb save$pairs
	cp pairs.*.rms  save$pairs
	cp pairs.*.plot save$pairs
end
