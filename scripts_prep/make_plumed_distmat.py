#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Bio.PDB
import sys
import glob
import numpy as np
import pandas as pd
import ast
skip_ev=3
fasta_file = open(sys.argv[1],mode='r')
PAE_cut=float(sys.argv[2])
Pr_cut=float(sys.argv[3])
print(PAE_cut,Pr_cut)
sequence = fasta_file.read().replace("\n","")
fasta_file.close()

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


#"MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGAGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPAKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRACDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSAGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGAGSASNAGSGSGFNGGFGSSMDSKSSGAGM"
file=open('resid_sel_all.dat', 'w')
#for i in range(1,len(sequence)+1,skip_ev):
for i in range(1,len(sequence)+1):
 #   for j in range(i+skip_ev,len(sequence),skip_ev):
    for j in range(i+skip_ev,len(sequence)):
        print(i,j)
        file.write(str(i)+" "+str(j)+"\n")
file.close()


# ### Read the AF csv contact map distances

# In[2]:


# relative path to search all text files
files = glob.glob("mean_af.csv")

#The contact map is in distance (Angstong)
print(files[0])
AF_contacts_file=files[0]
AF_resid_sel_file="resid_sel_all.dat"

AF_contacts=pd.read_csv(AF_contacts_file)
AF_contacts

exp=[]
labels=[]
lista=[]
residue_pairs=np.loadtxt(AF_resid_sel_file,usecols=(0,1))

####### Select residues that have a prob at p_d_i<0.03 at the last value

import json
import numpy as np
import matplotlib.pyplot as plt

js = glob.glob("*json")


f=open(js[0])
data=json.load(f)
l=0
PAE_array = np.empty(shape=(len(sequence),len(sequence)))
PAE_array.fill('NaN')

for i in data["predicted_aligned_error"]:
   print(i,len(i))
   for j in range(0,len(i)):
      #if (i[j]<5):
      PAE_array[l][j]=i[j]
      #print(i[j])
   print('l',l)
   l+=1

plt.matshow(PAE_array,cmap='jet',interpolation='none',vmin=0, vmax=40)
plt.colorbar()
plt.savefig('pae_m.png')
plt.close()

f.close()


pddf_pairs=np.load('prob.npy')
for i in range(len(residue_pairs)):

    r1=int(residue_pairs[i][0])
    r2=str(int(residue_pairs[i][1]))
    #This is to accoun that the trajectory starts residues from 1 while the AF contact map from 0.
    r1_m1=int(residue_pairs[i][0]-1)
    
    r1_label=str(r1)
    r2_label=r2
    ind1=int(r1_label)-1
    ind2=int(r2_label)-1

    if ( (PAE_array[ind1][ind2]+PAE_array[ind2][ind1])/2.0 <PAE_cut and (PAE_array[ind1][ind2]+PAE_array[ind2][ind1])/2.0 > 0 and pddf_pairs[ind1][ind2+1][-1]<Pr_cut):
        exp.append(AF_contacts.loc[r1_m1].at[r2])
        lista.append([int(r1_label),int(r2_label)])

########################
exp = np.array(exp)
labels=np.array(labels)

print(exp)
print (labels)
#len(exp)
print(len(lista))
print(lista)
with open('AF_contacts_constr.txt', 'w') as f:
    for line in range(len(exp)):
        f.write(str(exp[line])+"\n")
        


# ### AA af -> CG af

# In[ ]:

files = glob.glob("pdb_af.pdb")
print(files)
AF_pdb=files[0]
#"TDP_43_WtoA_7643e_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb"
name='input_af'

import mdtraj as md
traj=md.load(AF_pdb,top=AF_pdb)
selection_CA=traj.topology.select('name CA')

traj_CA=traj.atom_slice(selection_CA)

traj_CA.save(name+'.pdb')


topology = md.load(name+'.pdb').topology

table, bonds = topology.to_dataframe()

table['name'] = table['resName']
table['serial'] = table['resSeq']
table.loc[table.index[-1],'resName']='Z'
table.loc[table.index[0],'resName']='X'
table['name'] = table['resName']

print(table['resName'])

topologyCA = md.Topology.from_dataframe(table, bonds)
print(topologyCA)
traj = md.load_pdb(name+'.pdb',top=topologyCA)
traj.save(name+'.pdb')


# In[ ]:
######################
p = Bio.PDB.PDBParser()
#path = '/Volumes/BROTZAKIS/PROJECTS_LOCAL/AF_MI/Rela/prep' # your file path here
structure = p.get_structure('protein', AF_pdb)
plddt = [a.get_bfactor() for a in structure.get_atoms()]
aa=[a.get_name() for a in structure.get_atoms()]
res_ca=[]
for i in range(0,len(aa)):
        if (aa[i] == "CA"):
                res_ca.append(plddt[i])
plddt_cut=80
res_ca2=[n > plddt_cut for n in res_ca]
print(res_ca2)



print('key dis',ordered_domains)
print('key dis',disordered_domains)
######################
# In[ ]:



indexx=1
for i in ordered_domains.keys():
    print(ordered_domains[i],str(ordered_domains[i][0][0]),str(ordered_domains[i][0][1]))
    selection=traj.topology.select('resid '+str(ordered_domains[i][0][0])+' to '+str(ordered_domains[i][0][1]))
    print(selection)
    sel=traj.atom_slice(selection)
    sel.save('struct'+str(indexx)+'.pdb')
    indexx+=1


# In[7]:


f=open("plumed.dat","w")

f.write("MOLINFO MOLTYPE=protein STRUCTURE=input_af.pdb\n")
f.write("WHOLEMOLECULES ENTITY0=1-"+str(len(sequence))+"\n\n")
    
#Forward model
ordered_lista=[]
f.write("distance_rest_domains:  CONTACTMAP ...\n")
for i in range(0,len(lista)):
    for k in ordered_domains.keys():
        if ( lista[i][0] in [*range(ordered_domains[k][0][0],ordered_domains[k][0][1])] and lista[i][1] in [*range(ordered_domains[k][0][0],ordered_domains[k][0][1])]):
           ordered_lista.append(i)
ord_index=0
##########Make restrained contact map######
restr_array = np.empty(shape=(len(sequence),len(sequence)))
restr_array.fill('NaN')
for i in range(0,len(lista)):
    if (i not in ordered_lista):
        restr_array[lista[i][0]-1][lista[i][1]-1]=25

plt.matshow(restr_array,cmap='viridis_r',interpolation='none',vmin=0, vmax=25)
plt.colorbar()
plt.savefig('restr_distmap.png')
plt.close()



##################
for i in range(0,len(lista)):
    if (i not in ordered_lista):
        ord_index+=1
        f.write("ATOMS"+str(ord_index)+"="+str(lista[i][0])+","+str(lista[i][1])+"\n")
f.write("SWITCH={CUSTOM FUNC=x R_0=1}\n")
f.write("...\n\n")


#AF exp distance map definition. It has to be 1 to 1
f.write("af_dist_rest: CONSTANT VALUES=")
for i in range (0,len(exp)-1):
   if (i not in ordered_lista):
       f.write(str(exp[i]*0.1)+",")
       #f.write(str(exp[i])+",")

#This is converting to nm
if (len(exp)-1 not in ordered_lista):
   f.write(str(exp[len(exp)-1]*0.1)+" NODERIV\n\n")
f.write("af_dist2: CONSTANT VALUES=0 NODERIV\n")

#CV definition 
f.write("Rg: GYRATION TYPE=RADIUS ATOMS=1-"+str(len(sequence))+"\n")

indexx=1
cvnames=[]
for i in ordered_domains.keys():
    f.write("RMSD"+str(indexx)+": RMSD REFERENCE=struct"+str(indexx)+".pdb TYPE=OPTIMAL\n")
    cvnames.append("RMSD"+str(indexx))
    indexx+=1

indexx=1
RMSDs=[]
AT=[]
KAPPAS=[]
EXP=[]
EPS=[]
OFFSET=[]
for i in ordered_domains.keys():
    RMSDs.append("RMSD"+str(indexx))
    KAPPAS.append("100000")
    AT.append("0.02")
    EXP.append("2")
    EPS.append("1")
    OFFSET.append("0")
    indexx+=1

f.write("uwall: UPPER_WALLS ARG="+','.join(RMSDs)+" AT="+','.join(AT)+" KAPPA="+','.join(KAPPAS)+" EXP="+','.join(EXP)+" EPS="+','.join(EPS)+" OFFSET="+','.join(OFFSET)+"\n")


f.write("PRINT FILE=DISTANCE_MAP_REST ARG=distance_rest_domains.* STRIDE=200\n")
f.write("PRINT FILE=COLVAR ARG=__FILL__ STRIDE=200\n")
f.write("PBMETAD ...\n")
f.write("    LABEL=pb\n")
f.write("    ARG=__FILL__\n")
f.write("    SIGMA=__FILL__ \n")
f.write("    SIGMA_MIN=__FILL__ \n")
f.write("    SIGMA_MAX=__FILL__ \n")
f.write("    ADAPTIVE=DIFF\n")
f.write("    HEIGHT=__FILL__ \n")
f.write("    PACE=200\n")
f.write("    BIASFACTOR=__FILL__ \n")
f.write("    GRID_MIN=__FILL__ \n")
f.write("    GRID_MAX=__FILL__ \n")
f.write("    GRID_WSTRIDE=__FILL__ \n")
f.write("    WALKERS_MPI\n")
f.write("    TEMP=__FILL__ \n")
f.write("... PBMETAD\n")
f.write("METAINFERENCE ...\n")
f.write("    ARG=(distance_rest_domains.*),pb.bias REWEIGHT\n")
f.write("    PARARG=(af_dist_rest.*)\n")
f.write("    SIGMA_MEAN0=1\n")
f.write("    NOISETYPE=MGAUSS  OPTSIGMAMEAN=SEM AVERAGING=200\n")
f.write("    SIGMA0=10.0 SIGMA_MIN=0.0001 SIGMA_MAX=10.0 DSIGMA=0.1\n")
f.write("    MC_STEPS=10\n")
f.write("    MC_CHUNKSIZE=23\n")
f.write("    WRITE_STRIDE=10000\n")
f.write("    TEMP=__FILL__ \n")
f.write("    LABEL=af_mi_rest_domains\n")
f.write("... METAINFERENCE\n")
f.write("FLUSH STRIDE=200\n")
f.write("PRINT FILE=ENERGY ARG=pb.bias STRIDE=200\n")
f.write("PRINT ARG=af_mi_rest_domains.*   STRIDE=200 FILE=BAYES_rest_domains\n")
f.write("ENDPLUMED\n")

f.close()

