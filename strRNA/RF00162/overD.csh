# 1 = Ang. over target (17)

egrep '( A   9 | A  63 )'  save8/models.grem.pdb | awk -v a=$argv[1] '{if($2==9){x1=$7;y1=$8;y2=$9}else{x2=$7;y2=$8;y2=$9; x=x1-x2;y=y1-y2;z=z1-z2; d=sqrt(x*x+y*y+z*z); n++; if (d>17+a){m++}; print n,m}}' | tail -1
