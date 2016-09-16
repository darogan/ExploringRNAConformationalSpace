# 1 = RF00123, 2 = fraction gaps allowed (4)
# seqs.run aligns single sequences in align.seq
# alns.run aligns single/aligned files of sequences in align.lis

if ( -e seed.aln ) then
	echo seed alignment exists
else
	echo aligning seed sequences
	# rfam RF00123.afa.txt --> seq with gaps as "x"
	tcsh code/seed2seq.csh $argv[1]
	cat query.seq seed.seq > align.seq
	cat code/seqs.run | sed 's/50 500/50 2000/' > test.run
	./multas > multas.log
	grep seqs final.aln
	# delete SEED and convert fasta to seq
	sed -i '2 d' final.afa
	tcsh code/afa2seq.csh $argv[1]
	# final.afa written to seed.seq.cut (with gaps as "x")
	cat query.seq seed.cut.seq > final.seq
	echo aligning cut seed sequences
	./multas > multas.log
	cp final.aln seed.aln
	# check
	grep seqs final.aln
	grep '>' final.afa | wc -l
	echo top
	grep -v '>' final.afa | head
	echo bot
	grep -v '>' final.afa | tail
	cp final.afa seed.afa
endif
# combine seed and full sequences
#cat mulsel.seq seed.cut.seq | tr -d "x" > final.seq
#cp skip.seed.run skip.run
#./mulsel > mulsel.log
#tail -1 mulsel.log
#cp final.seq align.seq
cp mulsel.seq align.seq
echo aligning sequences with seed profile
cp code/align.lis .
if ( -e alns.run ) then
	cp alns.run test.run
else
	cp code/alns.run test.run
endif
./multas > multas.log
grep seqs final.aln
cp final.afa temp.afa
cp final.ali temp.ali
cp final.ali temp.tre
# remove short sequences (probably misaligned with the query)
tcsh code/cover.csh $argv[2]
# realign
echo realigning without rfam alignment
cp query.seq align.seq
awk '{if(match($1,">")){  if(match($0,"ALNS")){skip=1}else{skip=0} }; if(skip==0){print}}' final.afa > strip.afa
grep -v '>' strip.afa | tr "-" "x" | awk '{n++; print ">seq" n "\n" "seq\n" $1 "*"}' >> align.seq
# no gaps?
#grep -v '>' strip.afa | tr -d "-" | awk '{n++; print ">seq" n "\n" "seq\n" $1 "*"}' >> align.seq
# resort?
#grep -v '>' strip.afa | awk '{print $1, rand()}' | sort -n -k2 | awk '{print $1}' > strip.seq
#cat strip.seq | tr -d "-" | awk '{n++; print ">seq" n "\n" "seq\n" $1 "*"}' >> align.seq
if ( -e seqs.run ) then
	cp seqs.run test.run
else
	cp code/seqs.run test.run
endif
./multas > realign.log
grep seqs final.aln
cp align.seq full.seq
cp final.afa full.afa
cp final.ali full.ali
cp final.ali full.tre
# remove short sequences again
tcsh code/cover.csh $argv[2]
# check head and tail
grep '>' final.afa | wc -l
echo top
grep -v '>' final.afa | head
echo bot
grep -v '>' final.afa | tail
# check alignment
echo checking head and tail alignment
grep -v '>' final.afa | head | tr "-" "x" | awk '{n++; print ">seq" n "\n" "top\n" $1 "*"}' >  align.seq
grep -v '>' final.afa | tail | tr "-" "x" | awk '{n++; print ">seq" n "\n" "bot\n" $1 "*"}' >> align.seq
./multas > check.log
cat check.log | sed '1,/block 0 = 20 seqs/ d'
