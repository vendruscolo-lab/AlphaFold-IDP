n_files=$(ls ab_*.pdb|wc -l |awk '{print $1}')
for ((i=0;i<=$n_files;i++));do
name=ab_$i
../../src/pulchra -v $name".pdb"
file=$name"_em.pdb"
echo $name
if [ -f "$file" ]; then
        echo "$file exists"
else
#######  EM   ############
gmx_mpi pdb2gmx -f $name.rebuilt.pdb -o $name"_sys.pdb"  <<EOF 
6
1
EOF
gmx_mpi editconf -f $name"_sys.pdb" -o conf_box.gro -c -bt triclinic -d 1
gmx_mpi grompp -f  em.mdp  -c conf_box.gro -p topol.top -o em.tpr -maxwarn 28

gmx_mpi mdrun -deffnm em -v -g  -o em.trr
gmx_mpi trjconv -f em.gro -s em.gro -o $name"_em.pdb"<<EOF
0
EOF

rm topol.top *itp em.trr em.gro em.tpr em.log em.edr conf_box.gro $name"_sys.pdb"
###############
awk '{if (NR>5) print $0}'  $name"_em.pdb"   >>ab_bm_em.pdb
fi
done


gmx_mpi trjconv -f ab_bm_em.pdb -s ab_1_em.pdb -o abeta_C2.xtc <<EOF
0
EOF

cp ab_1_em.pdb topol_abeta_C2.pdb 
