#Make the plumed file
import mdtraj as md
import os
import sys
name='input_af'
pdb=name+'.pdb'
########### Conactenate ########
n_rep=int(sys.argv[1])
trajs = list()
skip=1
eq=1
os.system('rm input_af_*rebuilt.xtc input_af_bm_em.pdb input_af_bm_em_*.pdb *rebuilt.pdb segment*pdb *dcd')
for i in range(0,n_rep):
   try:
      print(f'loading {i}')
      traj = md.load_dcd(f'../output_{i}.dcd', top='input_af.pdb')

      traj.save('output_'+str(i)+'.xtc')
   except KeyboardInterrupt:
       raise
   except Exception as e:
       raise

