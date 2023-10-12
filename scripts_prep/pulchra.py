from mpi4py import MPI
import sys
import os
from glob import glob
import numpy as np

ntot=60000

k_file = open('k.dat',mode='r')
k = k_file.read().replace("\n","")
k_file.close()
k=int(k)
pdb_files = glob('./segment_'+str(k)+'_input_af_*.pdb')
print(pdb_files)
split=int(np.ceil(len(pdb_files)/5))
print ('SPLIT='+str(split))
start=0
end=int(np.ceil(len(pdb_files)))
print('end',end)
for j in range(start,len(pdb_files),2):
   print('In j loop', j,'lenpdb',len(pdb_files))
   name='segment_'+str(k)+'_input_af_'+str(j)
   os.system("/sharedscratch/fb516/pulchra304/src/pulchra -v {name}.pdb".format(name=name))
   os.system('echo "6 \n 7" | gmx pdb2gmx -f {name}.rebuilt.pdb -o {name}_sys.pdb -p topol_{name}.top -i {name}.itp'.format(name=name))
   os.system("rm topol_{name}.top {name}.itp em_{j}.trr em_{j}.tpr em_{j}.log em_{j}.edr conf_box_{j}.gro ".format(j=j,name=name))
   os.system('echo "0"|gmx trjconv -f  {name}_sys.pdb -s  {name}_sys.pdb -o {name}.rebuilt.xtc'.format(j=j,name=name))

