#jackhmmer -o jack.log -A jack.aln query.fa final.fa
nhmmer -E 0 -o jack.log -A jack.aln query.fa final.fa
cat jack.aln | grep USER | grep -v '=G[SRF]' | sed 's/.......................................... //' > jack.ali
cat jack.ali | awk -v n=1 '{print ">seq" n "\n" $1; n++}' | tr "." "-" > jack.faln
tcsh ~/TMpredict/codes/cutgap.csh jack 4
head jack.cut.faln
