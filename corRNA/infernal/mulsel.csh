# 1 = min kept (4000)

# set up mulsel run
@ n = `grep '>' final.seq | wc -l`
@ leng = `grep -v '[0-9\*]' pdb.seq | wc -c`
@ leng--
@ l  = $leng / 4
@ le = $leng - $l
@ ng = $leng + $l 
echo $n starting sequences
echo length = $leng
echo range = $le to $ng
cp code/skip99.run skip.run
eval "sed -i 's/50 500/$le $ng/' skip.run"
head -1 skip.run
echo running mulsel
code/mulsel > mulsel.log
tail -1 mulsel.log
# reduce with mulsel
foreach to ( 95 90 80 70 60 50 )
	echo reduce to $to
	cp code/skip$to.run skip.run
	code/mulsel > mulsel.log
	tail -1 mulsel.log
	cp final.seq final$to.seq
	@ n = `grep '>' final.seq | wc -l`
	if ( $n < $argv[1] ) break
end
@ half = $argv[1] / 2
if ( $n > $half ) then
	code/mulsel > mulsel.log
	tail -1 mulsel.log
	@ n = `grep '>' final.seq | wc -l`
endif
if ( $n > $half ) then
	code/mulsel > mulsel.log
	tail -1 mulsel.log
endif
