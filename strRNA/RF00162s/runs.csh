foreach move ( 0 2 4 6 8 10 )
	eval "cat code/stem.model | sed 's/MOVE/   $move/'" > stem.model
	echo >> runs.log
	echo move = $move >> runs.log
	tcsh code/pairs.csh 5 30 >> runs.log
	mkdir save$move
	cp models.*.pdb save$move
	cp pairs.*.rms  save$move
	cp pairs.*.plot save$move
end
