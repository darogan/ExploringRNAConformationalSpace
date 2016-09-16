# 1 = number of models

#@ n = `cat pairs.pred.dat | wc -l`

rm -r save*

@ len = `cat true.pdb | wc -l`
@ n = $len / 4 + 5
#@ n = ( 10 * $len ) / 35
echo using $n pairs
echo using $n pairs > runs.log
echo >> runs.log
@ pairs = $n + 1
while ( $pairs > 0 )
	@ pairs--
	echo pairs = $pairs
	eval "cat code/stem.model | sed 's/MOVE/   10/'" > stem.model
	echo >> runs.log
	echo pairs = $pairs >> runs.log
	tcsh code/pairs.csh 1 $pairs >> runs.log
	mkdir save$pairs
	cp models.*.pdb save$pairs
	cp pairs.*.rms  save$pairs
	cp pairs.*.plot save$pairs
end
