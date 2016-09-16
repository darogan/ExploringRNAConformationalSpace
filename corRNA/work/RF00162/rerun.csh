# 1 = min kept (4000)

cat home/infernal/RF00162/rfam.cut.afa | tr "U" "T" > rfam.afa

# set up mulsel run
cat rfam.afa | awk '{print $1; if(match($1,">")){print "rfam"}else{ print "*"}}' | sed 's/>seq1$/SEED:>seq1/' > final.seq
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
# reduce with mulsel
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
# compile reduced afa
grep seq final.seq | sed 's/>seq//' > mulsel.keep
head -2 rfam.afa > mulsel.afa
foreach seq ( `cat mulsel.keep | grep -v SEED` )
	@ at = $seq * 2
	head -$at rfam.afa | tail -2 >> mulsel.afa
end
@ n = `grep '>' mulsel.afa | wc -l`
echo $n sequences in mulsel.afa
