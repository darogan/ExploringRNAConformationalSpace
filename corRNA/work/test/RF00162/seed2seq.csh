cat $argv[1].afa.txt | sed 's/>/> /' | awk -v c=$argv[1] '{if($1==">"){n++; print "*\n>" c "=" n "\n" $2}else{print}}' | sed '1 d' | tr "-" "x" | sed 's/U/T/g' > seed.seq
echo "*" >> seed.seq
