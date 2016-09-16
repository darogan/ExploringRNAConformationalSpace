# 1 = number of models

rm -r save[0-9]
rm -r save1[0-9]

#@ n = `cat pairs.pred.dat | wc -l`
@ len = `cat true.pdb | wc -l`
@ n = $len / 4 + 5
echo using $n pairs
echo using $n pairs and $argv[1] models > runs.log
echo >> runs.log
foreach move ( 0 2 4 6 8 10 )
	@ moves = $move * 10000 + 1001
	echo Run length = $moves
	eval "cat code/config.run | sed 's/moves/ $moves/'" > config.run
	echo >> runs.log
	echo move = $move >> runs.log
	tcsh code/pairs.csh $argv[1] $n >> runs.log
	mkdir save$move
	cp models.*.pdb save$move
	cp pairs.*.rms  save$move
	cp pairs.*.plot save$move
end
