# 1 = min kept (4000)

cp query.seq final.seq
cat rfam.seq | sed 's/U/T/g' >> final.seq
@ n = `grep '>' final.seq | wc -l`
@ leng = `grep -v '[0-9\*]' query.seq | wc -c`
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
./mulsel > mulsel.log
tail -1 mulsel.log
#reduce
foreach to ( 95 90 80 70 60 50 )
	echo reduce to $to
	cp code/skip$to.run skip.run
	./mulsel > mulsel.log
	tail -1 mulsel.log
	cp final.seq final$to.seq
	@ n = `grep '>' final.seq | wc -l`
	if ( $n < $argv[1] ) break
end
@ half = $argv[1] / 2
if ( $n > $half ) then
	./mulsel > mulsel.log
	tail -1 mulsel.log
	@ n = `grep '>' final.seq | wc -l`
endif
if ( $n > $half ) then
	./mulsel > mulsel.log
	tail -1 mulsel.log
endif
#
echo mulsel selection copied to mulsel.seq
cat query.seq final.seq > mulsel.seq
# use addseed.csh now
exit

# align
echo aligning
cp query.seq mulsel.seq
grep -v 'SEED:' final.seq >> mulsel.seq
cp code/seqs.run .
./multas > multas.log
grep seqs final.aln
grep '>' final.afa | wc -l
echo top
grep -v '>' final.afa | head
echo bot
grep -v '>' final.afa | tail
