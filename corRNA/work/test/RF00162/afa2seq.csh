cat final.afa | sed 's/>/> /' | awk -v c=$argv[1] '{if($1==">"){n++; print "*\n>" c "=" n "\n" $2}else{print}}' | sed '1 d' | tr "-" "x" > seed.cut.seq
echo "*" >> seed.cut.seq
