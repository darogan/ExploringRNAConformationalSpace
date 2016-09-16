#jackhmmer -o jack.log -A jack.aln query.seq test.db
hmmbuild seed.hmm seed.afa
hmmalign seed.hmm final.fa > jack.aln
cat jack.aln | grep -v '=G[RSC]' | grep '>' | tr "." "-" | sed 's/.......................................... //' > jack.seq
cat jack.seq | awk -v n=1 '{print ">seq" n "\n" $1; n++}' > jack.faln
tcsh ~/TMpredict/codes/cutgap.csh jack 4
head jack.cut.faln
