# Supplementary Materials For Exploring RNA Conformational Space Under Sparse Distance Restraints #

## Contact Details ##
William R. Taylor & Russell S. Hamilton

## Repository Description ##

The two simulation programs are in simgen and simrna.
Links to the these will need to be changed to the local copies.

corRNA contains all the scripts used to set-up the sequence data.

strRNA contains the scripts to run the simgen simulations used to find the best perturbation level.
The top-level script is code/runs.csh which was run with a parameter 10 (making 10 models).
simRNA runs the same simulations for simrna.

strRNAfew runs the series of restraint removals.   In each application directory (eg RF00162)
executing the script code/runs.csh creats the full deletion series.    Similarly in simRNAfew


## Missing Files ##
Due to their large file size the following files are missing from the repository, but are publically available for download.

    corRNA/infernal/Rfam.cm (ftp://ftp.ebi.ac.uk/pub/databases/Rfam)
    corRNA/rfam9full/Rfam.full (ftp://ftp.ebi.ac.uk/pub/databases/Rfam)
    corRNA/save/seqs/pdb_seqres.txt (ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt)
    corRNA/seqs/pdb_seqres.txt (ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt)
