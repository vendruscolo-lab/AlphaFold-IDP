#!/bin/bash
nrep=$1
#<<BEG
rm COLVAR* HILLS* FULLBIAS
cp ../HILLS.* .
cp ../GRID* .
cp ../BAYES* .

rm segment*pdb segment*xtc
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

gmx trjcat -f nosolv_*.xtc -cat -o cat_trjcat.xtc -settime

#BEG
plumed driver --plumed plumed_analysis.dat --mf_xtc cat_trjcat.xtc
plumed --no-mpi driver --plumed reconstruct.dat --mf_xtc cat_trjcat.xtc --timestep 1

python resample.py
#Run the python script to make the fes plot
#num=1
#For other proteins the entries CV1,CV2,CV3 etc need to follow the COLVAR columns like:
#for i in $(echo CV1 CV2 CV3 etc);do

## For TDP-43 WtoA
#for i in $(echo Rg Rg1 Rg2 Rg3 Rg4 torsion1 torsion2 RMSD1 RMSD2 RMSD3);do
#python fes2.py --CV_col $num --CV_name $i
#num=$((num+1))
#echo $num
#done

python backmap.py $nrep
sh pulchra.sh
sh keepH.sh
