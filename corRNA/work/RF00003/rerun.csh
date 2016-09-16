# 1 = RF00123, 2 = min kept (4000)

cat home/infernal/$argv[1]/rfam.cut.afa | tr "U" "T" > rfam.afa

# set up mulsel run
cat rfam.afa | awk '{print $1; if(match($1,">")){print "rfam"}else{ print "*"}}' | sed 's/>seq1$/SEED:>seq1/' > final.seq
tcsh code/mulsel.csh $argv[2]
# compile reduced afa
grep seq final.seq | sed 's/>seq//' > mulsel.keep
head -2 rfam.afa > mulsel.afa
foreach seq ( `cat mulsel.keep | grep -v SEED` )
	@ at = $seq * 2
	head -$at rfam.afa | tail -2 >> mulsel.afa
end
@ n = `grep '>' mulsel.afa | wc -l`
echo $n sequences in mulsel.afa
