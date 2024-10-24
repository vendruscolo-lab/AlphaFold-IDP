#!/bin/bash
nrep=6
backmap=$1
#<<BEG
rm COLVAR* HILLS* FULLBIAS
cp ../HILLS.* .
cp ../GRID* .
cp ../BAYES* .

rm segment*pdb segment*xtc
python dcd2xtc.py
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
gmx trjcat -f nosolv_0.xtc  nosolv_1.xtc  nosolv_2.xtc  nosolv_3.xtc  nosolv_4.xtc  nosolv_5.xtc -cat -o cat_trjcat.xtc -settime <<EOF
0
c
c
c
c
c
EOF
#BEG
plumed driver --plumed plumed_analysis.dat --mf_xtc cat_trjcat.xtc
plumed --no-mpi driver --plumed reconstruct.dat --mf_xtc cat_trjcat.xtc --timestep 1

python resample.py
#Run the python script to make the fes plot
num=1
for i in $(echo CV1 CV2 CV3 etc);do
python fes2.py --CV_col $num --CV_name $i
num=$((num+1))
echo $num
done
python backmap.py
sh pulchra.sh
sh keepH.sh
