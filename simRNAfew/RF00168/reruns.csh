echo >> runs.log
echo rerun >> runs.log
foreach move ( 0 1 2 4 6 8 )
	eval "cat code/stem.model | sed 's/MOVE/   $move/'" > stem.model
	echo >> runs.log
	echo move = $move >> runs.log
	tcsh pairs.csh 10 60 >> runs.log
	#mkdir save$move
	cp models.*.pdb save$move
	cp pairs.*.rms  save$move
	cp pairs.*.plot save$move
end
