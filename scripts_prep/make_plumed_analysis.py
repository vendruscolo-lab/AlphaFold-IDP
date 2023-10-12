#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Bio.PDB
import sys
import glob
import numpy as np
import pandas as pd
import ast


df = pd.read_csv('ordered_domains.csv')
kkeys=[]
lis=[]
# Convert the string column back to a list column
for i in df.keys():
    hi = df[i].apply(ast.literal_eval)
    kkeys.append(i)
    lis.append(hi[0])

print(lis)
my_df  = pd.DataFrame(columns = kkeys)
my_df.loc[len(my_df)] = lis
print(my_df)
ordered_domains=my_df.to_dict('list')


df = pd.read_csv('disordered_domains.csv')
kkeys=[]
lis=[]
# Convert the string column back to a list column
for i in df.keys():
    hi = df[i].apply(ast.literal_eval)
    kkeys.append(i)
    lis.append(hi[0])

print(lis)
my_df  = pd.DataFrame(columns = kkeys)
my_df.loc[len(my_df)] = lis
print(my_df)
disordered_domains=my_df.to_dict('list')


fasta_file = open(sys.argv[1],mode='r')
sequence = fasta_file.read().replace("\n","")
fasta_file.close()

print('key dis',ordered_domains)
print('key dis',disordered_domains)




f=open("plumed_analysis.dat","w")
r=open("reconstruct.dat","w")
f.write("RESTART\n")
f.write("MOLINFO MOLTYPE=protein STRUCTURE=input_af.pdb\n")
r.write("MOLINFO MOLTYPE=protein STRUCTURE=input_af.pdb\n")
f.write("WHOLEMOLECULES ENTITY0=1-"+str(len(sequence))+"\n\n")
r.write("WHOLEMOLECULES ENTITY0=1-"+str(len(sequence))+"\n\n")
r.write("DUMPATOMS STRIDE=1 FILE=conf_0_recon.gro ATOMS=1-"+str(len(sequence))+"\n")

r.close()

cvnames=[]
#CV definition 
f.write("Rg: GYRATION TYPE=RADIUS ATOMS=1-"+str(len(sequence))+"\n")
cvnames.append("Rg")
indexx=1
for i in ordered_domains.keys():
    f.write("RMSD"+str(indexx)+": RMSD REFERENCE=struct"+str(indexx)+".pdb TYPE=OPTIMAL\n")
    cvnames.append("RMSD"+str(indexx))
    indexx+=1
indexx=1


indexx=1

f.write("t1: CENTER ATOMS=3-40\n")
f.write("t2: CENTER ATOMS=41-79\n")
f.write("t3: CENTER ATOMS=104-140\n")
f.write("t4: CENTER ATOMS=141-178\n")

f.write("torsion1: TORSION ATOMS=t1,t2,t3,t4\n")

f.write("tt1: CENTER ATOMS=104-140\n")
f.write("tt2: CENTER ATOMS=141-178\n")
f.write("tt3: CENTER ATOMS=191-225\n")
f.write("tt4: CENTER ATOMS=226-260\n")

f.write("torsion2: TORSION ATOMS=tt1,tt2,tt3,tt4\n")


f.write("# PBMetaD\n")
f.write("PBMETAD ...\n")
f.write("    LABEL=pb\n")
f.write("    ARG=torsion1,torsion2\n")
f.write("    SIGMA=1000\n")
f.write("    SIGMA_MIN=0.05,0.05\n")
f.write("    SIGMA_MAX=0.1,0.1\n")
f.write("    ADAPTIVE=DIFF\n")
f.write("    HEIGHT=0.5\n")
f.write("    PACE=20000000\n")
f.write("    BIASFACTOR=35\n")
f.write("    GRID_MIN=-pi,-pi\n")
f.write("    GRID_MAX=pi,pi\n")
f.write("    GRID_WSTRIDE=5000\n")
f.write("    WALKERS_MPI\n")
f.write("    TEMP=298\n")
f.write("... PBMETAD\n")


cv=','.join(cvnames)
f.write("PRINT ARG=pb.bias FILE=FULLBIAS STRIDE=1\n")
f.write("PRINT FILE=COLVAR ARG="+str(cv)+" STRIDE=1\n")
f.write("ENDPLUMED\n")

f.close()

