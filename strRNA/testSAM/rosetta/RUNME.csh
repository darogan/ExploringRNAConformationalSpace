# 1 = file(.pdb)

# convert to "DNA"
cat ../2gis.pdb | grep '^ATOM ' | sed 's/ U A / T A /' > 2gis.pdb

# add Hs
main/fixbb.static.linuxgccrelease -database main/minirosetta_database -packing:ex1 -packing:use_input_sc -in:file:s $argv[1].pdb -resfile empty.resfile
mv $argv[1]_0001.pdb $argv[1].H.pdb

# count H-bonde
cat  $argv[1].H.pdb | grep '^ATOM ' | grep -v ' [CP]' > bonds.dat
grep    H bonds.dat | cut -c 23-99 > hydro.dat
grep -v H bonds.dat | cut -c 23-99 > heavy.dat

