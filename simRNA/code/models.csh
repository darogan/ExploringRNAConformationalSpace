# 1 = runs
 
set pdb = `head -1 full.pdb`
if ( $pdb[5] != "A" ) then
	echo Chain must be A
	exit
endif
if ( $pdb[6] != "1" ) then
	echo Chain must start at 1
	exit
endif
awk '{printf("SLOPE A/%d/MB A/%d/MB 5.0 10.0 1.0\nWELL A/%d/MB A/%d/MB 5.0 10.0 5.0\n", $1,$2,$1,$2)}' pairs.dat > pairs.con

if ( -e models.pdb ) rm models.pdb
@ n = 0
while ( $n < $argv[1] )
	@ n++
	if ( -e pred.seq ) then
		# predict
		code/sim -s pred.seq -o models -c config.run -S bracket.dat -r pairs.con >& run.log
	else
		# denature
		code/sim -p full.pdb -o models -c config.run -S bracket.dat -r pairs.con >& run.log
	endif
	rm models-*.pdb
	set run = `head -1 config.run`
	@ get = 1 + $run[2] / 1000
	code/sim2pdb full.pdb models.trafl $get >& /dev/null
	tcsh code/super.csh 1
	grep ' A ' super.pdb >> models.pdb
	echo TER >> models.pdb
end
# add true.pdb and view 
grep ' B ' super.pdb >> models.pdb
# rasmol -script code/models.ras models.pdb &
