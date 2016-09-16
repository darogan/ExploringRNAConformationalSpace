# 1 = fraction of length

echo removing sequences with more than 1/$argv[1] gaps
grep -v '>' final.afa > final.tmp
rm final.afa
@ len = `head -1 final.tmp | wc -c`
@ max = $len / $argv[1]
foreach seq ( `cat final.tmp` )
        @ left = `echo $seq | tr -d "-" | wc -c`
        @ gaps = $len - $left
        if ( $gaps > $max ) continue
        echo $seq | awk '{n++; print ">seq" n "\n" $1}' >> final.afa
end
