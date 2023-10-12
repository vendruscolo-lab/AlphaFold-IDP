#!/usr/bin/env python
# coding: utf-8



#Make the plumed file
import mdtraj as md
import os
name='input_af'
pdb=name+'.pdb'
########### Conactenate ########
n_rep=6
trajs = list()
skip=2
eq=10


#############
xtc='r0/conf-protein.xtc'
topology = md.load(pdb).topology
traj = md.load(xtc, top=pdb)

table, bonds = topology.to_dataframe()
print('CALVADOS topology',table.head())

#########
fasta_file = open('sequence.dat',mode='r')
fasta = fasta_file.read().replace("\n","")
fasta_file.close()
print(fasta[0])

aminoacids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
aminos_sin=['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

print(aminoacids[aminos_sin.index(fasta[0])])
print(aminoacids[aminos_sin.index(fasta[-1])])
table.loc[table.index[-1],'resName']=aminoacids[aminos_sin.index(fasta[-1])]
table.loc[table.index[0],'resName']=aminoacids[aminos_sin.index(fasta[0])]

table.loc[table.index[-1],'name']=aminoacids[aminos_sin.index(fasta[-1])]
table.loc[table.index[0],'name']=aminoacids[aminos_sin.index(fasta[0])]
table['serial'] = table['resSeq']
########

for AA in aminoacids:
    table['name'] = table['name'].replace([AA], 'CA')


print('new topology with CA, cg naming', table.head())


topologyCA = md.Topology.from_dataframe(table, bonds)


print(topologyCA)


traj = md.load_pdb(pdb,top=topologyCA)
traj.save(name+'cg.pdb')




#Trajectory
for  j in range(5,5+1):
   name='input_af'
   pdb=name+'.pdb'
   dcd='r0/conf-protein.xtc'
   traj_name=name+'_'+str(j)+'dcd.pdb'
   topology = md.load(pdb).topology
   traj = md.load(dcd, top=pdb)
   traj.save(traj_name)
   pdb=traj_name
   table, bonds = topology.to_dataframe()
   print('CALVADOS topology',table.head())
   #########
   fasta_file = open('sequence.dat',mode='r')
   fasta = fasta_file.read().replace("\n","")
   fasta_file.close()
   print(fasta[0])
 
   aminoacids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
   aminos_sin=['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

   print(aminoacids[aminos_sin.index(fasta[0])])
   print(aminoacids[aminos_sin.index(fasta[-1])])
   table.loc[table.index[-1],'resName']=aminoacids[aminos_sin.index(fasta[-1])]
   table.loc[table.index[0],'resName']=aminoacids[aminos_sin.index(fasta[0])]

   table.loc[table.index[-1],'name']=aminoacids[aminos_sin.index(fasta[-1])]
   table.loc[table.index[0],'name']=aminoacids[aminos_sin.index(fasta[0])]
   table['serial'] = table['resSeq']
   ########
   for AA in aminoacids:
      table['name'] = table['name'].replace([AA], 'CA')
   print('new topology with CA, cg naming', table.head())
   topologyCA = md.Topology.from_dataframe(table, bonds)
   traj = md.load_pdb(pdb,top=topologyCA)
   traj.save(name+'_'+str(j)+'dcdcg.pdb')
   print('laste'+str(j))
   for i in range(0,len(traj)):
      traj[i].save('segment_'+str(j)+'_'+name+'_'+str(i)+'.pdb')
