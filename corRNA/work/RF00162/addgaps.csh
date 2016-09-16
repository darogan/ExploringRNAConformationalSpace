rm mulsel.gap.afa
foreach line ( `cat mulsel.fa` )
	set code = `echo $line | sed 's/>RFAM/R /'`
	if ( $code[1] == 'R' ) then
		# enter code and gapped seq
		echo $line >> mulsel.gap.afa
		eval "grep '^$code[2]	' home/rfam9full/RF00162/RF00162.cut.xgaps" | awk '{print $2}' >> mulsel.gap.afa
		echo $line >> mulsel.gap.afa
	else
		echo $line >> mulsel.gap.afa
	endif
end 
