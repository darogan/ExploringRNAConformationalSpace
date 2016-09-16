cat blast162.txt | awk  '{n++; if($1=="Query"){ n=0 } print n, $1, $3}' > blast.tmp1
sort -k1 -n blast.tmp1 | awk '{print $2, $3}' > blast.tmp2
cat blast.tmp2 | awk '{if($1==last){print $2}else{print "#>" $1 "#" "\n" $2; last = $1}}' | tr -d "\n" | tr "#" "\n" > blast162.afa
more blast162.afa
