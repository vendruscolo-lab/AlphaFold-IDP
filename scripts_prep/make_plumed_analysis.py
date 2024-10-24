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

f.write("# PBMetaD\n")
f.write("PBMETAD ...\n")
f.write("    LABEL=pb\n")
f.write("    ARG=__FILL__\n")
f.write("    SIGMA=__FILL__\n")
f.write("    SIGMA_MIN=__FILL__\n")
f.write("    SIGMA_MAX=__FILL__\n")
f.write("    ADAPTIVE=DIFF\n")
f.write("    HEIGHT=__FILL__\n")
f.write("    PACE=20000000\n")
f.write("    BIASFACTOR=__FILL__\n")
f.write("    GRID_MIN=__FILL__\n")
f.write("    GRID_MAX=__FILL__\n")
f.write("    GRID_WSTRIDE=__FILL__\n")
f.write("    WALKERS_MPI\n")
f.write("    TEMP=__FILL__\n")
f.write("... PBMETAD\n")


cv=','.join(cvnames)
f.write("PRINT ARG=pb.bias FILE=FULLBIAS STRIDE=1\n")
f.write("PRINT FILE=COLVAR ARG=__FILL__ STRIDE=1\n")
f.write("ENDPLUMED\n")

f.close()

