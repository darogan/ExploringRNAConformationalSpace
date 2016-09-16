# 1 = min kept (4000), 2 = range +/- fraction, default (4)

# set up mulsel run
@ n = `grep '>' final.seq | wc -l`
@ leng = `grep -v '[0-9\*]' pdb.seq | wc -c`
@ leng--
@ l  = $leng / $argv[2]
@ le = $leng - $l
@ ng = $leng + $l 
echo $n starting sequences
echo length = $leng
echo range = $le to $ng
cp code/skip99.run skip.run
eval "sed -i 's/50 500/$le $ng/' skip.run"
head -1 skip.run
echo running mulsel
cp final.seq start.seq
code/mulsel > mulsel.log
tail -1 mulsel.log
cp final.seq final99.seq
# reduce with mulsel
foreach to ( 95 90 85 80 70 60 50 )
	@ n = `grep '>' final.seq | wc -l`
	if ( $n < $argv[1] ) break
	echo reduce to $to
	cp code/skip$to.run skip.run
	code/mulsel > mulsel.log
	tail -1 mulsel.log
	cp final.seq final$to.seq
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
