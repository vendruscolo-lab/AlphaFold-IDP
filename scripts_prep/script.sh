#!/bin/bash
nrep=$1
rm COLVAR* HILLS* FULLBIAS
cp ../HILLS.* .
cp ../GRID* .
cp ../BAYES* .

rm segment*pdb segment*xtc
#Convert dcds to xtc
python dcd2xtc.py $nrep
for ((c=0;c<$nrep;c++))
do

for part in $(ls output_"$c".xtc)
do
echo $part
gmx trjconv -f $part -s ../input_af.pdb -o nosolv_$c.xtc <<EOF
0
EOF
done
done

#concatenate xtcs to a single xtc
gmx trjcat -f nosolv_*.xtc -cat -o cat_trjcat.xtc -settime

plumed driver --plumed plumed_analysis.dat --mf_xtc cat_trjcat.xtc
plumed --no-mpi driver --plumed reconstruct.dat --mf_xtc cat_trjcat.xtc --timestep 1

#Sample the structural ensemble by weights.
python resample.py


#Backmap from coarse-grained to atomistic
python backmap.py $nrep
sh pulchra.sh
sh keepH.sh
