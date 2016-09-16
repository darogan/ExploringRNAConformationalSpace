# 1 = number of models

#@ n = `cat pairs.pred.dat | wc -l`
@ len = `cat true.pdb | wc -l`
@ n = $len / 4 + 5
echo using $n pairs
#echo using $n pairs and $argv[1] models > runs.log
#echo >> runs.log
foreach move ( 6 8 10 )
	eval "cat code/stem.model | sed 's/MOVE/   $move/'" > stem.model
	echo >> runs.log
	echo move = $move >> runs.log
	tcsh code/pairs.csh $argv[1] $n >> runs.log
	mkdir save$move
	cp models.*.pdb save$move
	cp pairs.*.rms  save$move
	cp pairs.*.plot save$move
end
