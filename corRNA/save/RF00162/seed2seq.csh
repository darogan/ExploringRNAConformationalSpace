cat RF00162.afa.txt | sed 's/>/> /' | awk '{if($1==">"){print "*\nS-" $1 $2 "\n" $2}else{print}}' | sed '1 d' | tr "-" "x" | sed 's/U/T/g' > seed.seq
echo "*" >> seed.seq
