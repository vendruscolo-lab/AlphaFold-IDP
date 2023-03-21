#%%bash 
##Make constraint input. 
#rm /Volumes/BROTZAKIS/PROJECTS_LOCAL/ALPHA_IDP/KULL_Centre_papers_main_2022_rh_fwd_model_pesce_et_al/ALL_DISTANCES_AF/distmat_with_distributions/Tau_distmat/resid_sel_all.dat
#for ((c=1;c<=441-1;c+=4))do
# for ((j=$c+4;j<=441;j+=4))do
#  echo $c $j >>  /Volumes/BROTZAKIS/PROJECTS_LOCAL/ALPHA_IDP/KULL_Centre_papers_main_2022_rh_fwd_model_pesce_et_al/ALL_DISTANCES_AF/distmat_with_distributions/Tau_distmat/resid_sel_all.dat
#done 
#done
#rm  /Volumes/BROTZAKIS/PROJECTS_LOCAL/ALPHA_IDP/KULL_Centre_papers_main_2022_rh_fwd_model_pesce_et_al/ALL_DISTANCES_AF/distmat_with_distributions/Tau_distmat/resid_sel_all_ev1.dat
#for ((c=1;c<=441-1;c+=2))do
# for ((j=$c+2;j<=441;j+=2))do
#  echo $c $j >>  /Volumes/BROTZAKIS/PROJECTS_LOCAL/ALPHA_IDP/KULL_Centre_papers_main_2022_rh_fwd_model_pesce_et_al/ALL_DISTANCES_AF/distmat_with_distributions/Tau_distmat/resid_sel_all_ev1.dat
#done 
#done
#######

import sys,argparse,os
import pandas as pd
import numpy as np
import subprocess
import glob
import matplotlib.pyplot as plt
import matplotlib
import mdtraj as md
import argparse, sys
import statistics

parser = argparse.ArgumentParser()
parser.add_argument("--master_folder", help="path where csv and output_files files lie/wil lie",required=True,type=str)
parser.add_argument("--pdb", help="pdb path",required=True,type=str)
parser.add_argument("--xtc", help="xtc path",required=True,type=str)
parser.add_argument("--sys_n", help="system name",required=True,type=str)
parser.add_argument("--SAXS_PPDF_file", help="SAXS PDDF file out file",required=True,type=str)
parser.add_argument("--ev_CB", help="every CB",required=True,type=int)
parser.add_argument("--nresid", help="number of residues",required=True,type=int)
parser.add_argument("--SAXS_bw", help="times 10 for SAXS?",required=True,type=int)

args = parser.parse_args()
data_f =args.master_folder+"/"
pdb =args.pdb
xtc =args.xtc
system_name =args.sys_n
SAXS_PPDF_file =args.SAXS_PPDF_file
ev_CB =args.ev_CB
nresid=args.nresid
SAXS_bw=args.SAXS_bw


def weighted_median(data, weights):
    """
    Args:
      data (list or numpy.array): data
      weights (list or numpy.array): weights
    """
    data, weights = np.array(data).squeeze(), np.array(weights).squeeze()
    s_data, s_weights = map(np.array, zip(*sorted(zip(data, weights))))
    midpoint = 0.5 * sum(s_weights)
    if any(weights > midpoint):
        w_median = (data[weights == np.max(weights)])[0]
    else:
        cs_weights = np.cumsum(s_weights)
        idx = np.where(cs_weights <= midpoint)[0][-1]
        if cs_weights[idx] == midpoint:
            w_median = np.mean(s_data[idx:idx+2])
        else:
            w_median = s_data[idx+1]
    return w_median

file=open(data_f+'resid_sel_all.dat', 'w')

for i in range(1,nresid,1):
    for j in range(i+ev_CB,nresid+1,1):
        file.write(str(i)+" "+str(j)+"\n")
file.close()        
#The contact map is in distance (Angstong)
#data_f="/Volumes/BROTZAKIS/PROJECTS_LOCAL/ALPHA_IDP/KULL_Centre_papers_main_2022_rh_fwd_model_pesce_et_al/ALL_DISTANCES_AF/distmat_with_distributions/Tau_distmat/"
AF_contacts_file=data_f+"alphafold2_ptm_model_3_seed_000_mean.csv"
AF_resid_sel_file=data_f+"resid_sel_all.dat"
#AF_resid_sel_file=data_f+"resid_sel_all_ev1.dat"
AF_contacts_stdev_file=data_f+"alphafold2_ptm_model_3_seed_000_std.csv"
AF_contacts=pd.read_csv(AF_contacts_file)
AF_contacts_stdev=pd.read_csv(AF_contacts_stdev_file)
AF_contacts
AF_contacts_stdev

exp=[]
labels=[]
lista=[]
residue_pairs=np.loadtxt(AF_resid_sel_file,usecols=(0,1))
for i in range(len(residue_pairs)):

    r1=int(residue_pairs[i][0])
    r2=str(int(residue_pairs[i][1]))
    #This is to accoun that the trajectory starts residues from 1 while the AF contact map from 0.
    r1_m1=int(residue_pairs[i][0]-1)
    
    r1_label=str(r1)
    r2_label=r2
    
    #The first index of loc needs to be integer and minus1, since it starts from 0. index 2 needs to be str and reads properly
    exp.append([AF_contacts.loc[r1_m1].at[r2],AF_contacts_stdev.loc[r1_m1].at[r2]])
    #exp.append([AF_contacts.loc[r1_p1].at[r2_p1],2])
    labels.append('resid_'+r1_label+'-'+r2_label)

    lista.append([int(r1_label),int(r2_label)])
    #lista.append([r1_p1,r2_p1])

exp = np.array(exp)

labels=np.array(labels)
labels
print(exp)
print (labels)
len(exp)
print(len(lista))
print(lista)
with open(data_f+'AF_contacts_constr.txt', 'w') as f:
    f.write("# DATA=JCOUPLINGS\n")
    for line in range(len(exp)):
        f.write(str(labels[line])+" "+str(exp[line,0])+" "+str(exp[line,1])+"\n")
        
######
#Make the plumed file
#pdb=data_f+'../../../trajectories/Tau/0_em.pdb'
#xtc=data_f+'../../../trajectories/Tau/trjout_bm_em.xtc'

topology = md.load(pdb).topology
traj = md.load(xtc, top=pdb, stride=50)

#traj=traj1[::50]
#traj1=[]
lista_CB=[]
for i in lista:
    if (topology.select("residue "+str(i[0])+" and name CA and (resname  GLY PRO)")):
        res1=topology.select("residue "+str(i[0])+" and name CA and (resname  GLY PRO)")
        #print(i[0],res1)
    if (topology.select("residue "+str(i[1])+" and name CA and (resname  GLY PRO)")):
        res2=topology.select("residue "+str(i[1])+" and name CA and (resname  GLY PRO)")
        #print(i[1],res2)
    if (topology.select("residue "+str(i[0])+" and name CB")):
        res1=topology.select('name CB and residue '+str(i[0]))
        #print(i[0],res1)
    if (topology.select("residue "+str(i[1])+" and name CB")):
        res2=topology.select('name CB and residue '+str(i[1]))
        #print(i[1],res2)
    lista_CB.append([int(res1),int(res2)])

print('CB-atomlist')
for j in lista_CB:
    print(j[0],j[1])
print(len(lista_CB))
print(lista_CB)

for i in lista:
    print(i)


distance_nm=md.compute_distances(traj,lista_CB)


#This is converting to Angstorng
calc=distance_nm[:] * 10
print(calc)

with open(data_f+'AF_contacts_calc.txt', 'w') as f:
    for line in range(len(calc[:,0])):

        f.write(str(line)+" ")
        for j in range(len(lista)): 
            f.write(str(calc[line,j])+" ")
        
        f.write("\n")        
        
#print(calc[:,0])
################

#print(len(calc[:,1]))
#print(calc[:,1])
file=open(data_f+'CG_distances.dat', 'w')
file.write("\n")
for i in range(0,calc.shape[1]):
    file.write(str(statistics.median(calc[:,i]))+"\n")
file.close()
file_name=data_f+'AF_CG.dat'
lis=subprocess.call(['paste'] + glob.glob(data_f+'AF_contacts_constr.txt') + [data_f+'CG_distances.dat'] , stdout=open(file_name, 'w') )
MD=np.loadtxt(file_name,usecols=(3))
AF=np.loadtxt(file_name,usecols=(1))
#####

plt.scatter(MD, AF, marker='+',lw=1)
#plt.xlim(0,20)
plt.xlim(5,21.84)
plt.ylim(5,21.84)
plt.xlabel("CG distance (A)",size=15) 
plt.ylabel("AF distance (A)",size=15) 
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#plt.axline((1, 1), slope=1,ls='--',color='black')
#plt.xlim(5,20)
#plt.ylim(5,20)
plt.savefig(data_f+'CG_AF.pdf',bbox_inches='tight')
plt.clf()
plt.close()
###
file_name=data_f+'AF_CG.dat'
lis=subprocess.call(['paste'] + glob.glob(data_f+'AF_contacts_constr.txt') + [data_f+'CG_distances.dat'] , stdout=open(file_name, 'w') )
MD=np.loadtxt(file_name,usecols=(3))
AF=np.loadtxt(file_name,usecols=(1))
#####

plt.scatter(MD, AF, marker='+',lw=1)
#plt.xlim(0,20)
#plt.xlim(5,21.84)
#plt.ylim(5,21.84)
plt.xlabel("CG distance (A)",size=15)
plt.ylabel("AF distance (A)",size=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#plt.axline((1, 1), slope=1,ls='--',color='black')
#plt.xlim(5,20)
#plt.ylim(5,20)
plt.savefig(data_f+'CG_AF_long.pdf',bbox_inches='tight')
plt.clf()
plt.close()


rg_A=md.compute_rg(traj)*10



#Experimental PDDF is in nm.
pddf_SAXS_f=SAXS_PPDF_file
data = np.loadtxt(pddf_SAXS_f)
#*10 converts to anstrong
ppdf_SAXS_bin=data[:,0]*SAXS_bw
ppdf_SAXS_prob=data[:,1]
SAXS_err=data[:,2]

#list1 = [list2[i//n] for i in range(len(list1))]

plt.title(system_name)
plt.plot(ppdf_SAXS_bin, ppdf_SAXS_prob/sum(ppdf_SAXS_prob),c='black',label="SAXS",lw=1)
plt.scatter(ppdf_SAXS_bin, ppdf_SAXS_prob/sum(ppdf_SAXS_prob), marker='+',lw=1,c='black')
pddf_pairs=np.load(data_f+'alphafold2_ptm_model_3_seed_000_prob_distributions.npy')


#Ignoring last bin
ppdf_AF_prob=[0.0]*(len(pddf_pairs[0][1])-1)
ppdf_AF_bin=[]
AF_bin_size=0.3125
for i in range(0,len(pddf_pairs[0][1])-1):
    ppdf_AF_bin.append(i*AF_bin_size)

print(len(ppdf_AF_prob))

print(len(pddf_pairs[0][1]))
for i in range(0,len(pddf_pairs[0])-1):
    for j in range(i+ev_CB,len(pddf_pairs[0])):
        for k in range(len(pddf_pairs[0][1])-1):
            #print(i, j,k, pddf_pairs[i][j][k])
            ppdf_AF_prob[k]+=pddf_pairs[i][j][k]
plt.xlim(0,20)

#Integral approach: 
norm_ppdf_AF_prob=ppdf_AF_prob/sum(ppdf_AF_prob)
norm_ppdf_SAXS_prob=ppdf_SAXS_prob/sum(ppdf_SAXS_prob)
sum_AF=sum(norm_ppdf_AF_prob[0:len(ppdf_AF_bin)])
binwidth_SAXS=(ppdf_SAXS_bin[1]-ppdf_SAXS_bin[0])
bin_SAXS_cut=int(20/binwidth_SAXS)
sum_SAXS=sum(norm_ppdf_SAXS_prob[0:bin_SAXS_cut])
norm=(sum_AF*AF_bin_size)/(sum_SAXS*binwidth_SAXS)

plt.plot(ppdf_AF_bin,norm_ppdf_AF_prob/norm,c='b',label='AF')
plt.scatter(ppdf_AF_bin, norm_ppdf_AF_prob/norm, marker='+',lw=1,c='b')

plt.xlabel("r ($\AA$)",fontsize=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.ylabel('P',fontsize=15)
plt.legend()
plt.savefig(data_f+'PDDF.pdf',bbox_inches='tight')
plt.clf()
plt.close()


resid=list(np.loadtxt(file_name,usecols=(0),dtype=str))
file=open(data_f+"resid_sel_all_sample.dat", 'w')
#Sample random from AF PDDF
print(sum(norm_ppdf_AF_prob))
new_random=np.random.choice(ppdf_AF_bin,100,p=norm_ppdf_AF_prob)
print(new_random)
res_lt1=[]
past=[]
list(set(AF).intersection(new_random))
for j in range(0,len(new_random)):
    for i in range(0,len(AF)):
        #if (np.abs(new_random[j]-AF[i])<0.3 and MD[i]<20):
        #if (np.abs(new_random[j]-AF[i])<0.3 and MD[i]<20 and AF[i]<18):
        if (np.abs(new_random[j]-AF[i])<0.3 and AF[i]<16 and np.abs(MD[i]-AF[i])<3.1):
            print (AF[i],new_random[j],j,i,MD[i])
            if (i not in past):
                res_lt1.append(str(resid[i].split("-")[0].split("_")[1])+" "+str(resid[i].split("-")[1]))
                past.append(i)
                break

    #print(AF[i])
#for i in range(0,len(AF)):
#    if (AF[i]<17 and MD[i]<17):
#    #if (MD[i]<20):
#        res_lt1.append(str(resid[i].split("-")[0].split("_")[1])+" "+str(resid[i].split("-")[1]))
        

print(len(res_lt1))
for i in range(0,len(res_lt1)):
    file.write(str(res_lt1[i])+"\n")

file.close()


#The contact map is in distance (Angstong)
#data_f="/Volumes/BROTZAKIS/PROJECTS_LOCAL/ALPHA_IDP/KULL_Centre_papers_main_2022_rh_fwd_model_pesce_et_al/ALL_DISTANCES_AF/distmat_with_distributions/Tau_distmat/"
AF_contacts_file=data_f+"alphafold2_ptm_model_3_seed_000_mean.csv"
AF_resid_sel_file=data_f+"resid_sel_all_sample.dat"
AF_contacts_stdev_file=data_f+"alphafold2_ptm_model_3_seed_000_std.csv"

AF_contacts=pd.read_csv(AF_contacts_file)
AF_contacts_stdev=pd.read_csv(AF_contacts_stdev_file)
AF_contacts
AF_contacts_stdev

exp=[]
labels=[]
lista=[]
residue_pairs=np.loadtxt(AF_resid_sel_file,usecols=(0,1))
for i in range(len(residue_pairs)):

    r1=int(residue_pairs[i][0])
    r2=str(int(residue_pairs[i][1]))
    #This is to accoun that the trajectory starts residues from 1 while the AF contact map from 0.
    r1_m1=int(residue_pairs[i][0]-1)
    
    r1_label=str(r1)
    r2_label=r2
    
    #The first index of loc needs to be integer and minus1, since it starts from 0. index 2 needs to be str and reads properly
    exp.append([AF_contacts.loc[r1_m1].at[r2],AF_contacts_stdev.loc[r1_m1].at[r2]])
    #exp.append([AF_contacts.loc[r1_p1].at[r2_p1],2])
    labels.append('resid_'+r1_label+'-'+r2_label)

    lista.append([int(r1_label),int(r2_label)])
    #lista.append([r1_p1,r2_p1])

exp = np.array(exp)

labels=np.array(labels)
labels
print(exp)
print (labels)
len(exp)
print(len(lista))
print(lista)
with open(data_f+'AF_contacts_constr.txt', 'w') as f:
    f.write("# DATA=JCOUPLINGS\n")
    for line in range(len(exp)):
        f.write(str(labels[line])+" "+str(exp[line,0])+" "+str(exp[line,1])+"\n")


lista_CB=[]      
for i in lista:
    if (topology.select("residue "+str(i[0])+" and name CA and (resname  GLY PRO)")):
        res1=topology.select("residue "+str(i[0])+" and name CA and (resname  GLY PRO)")
        #print(i[0],res1)
    if (topology.select("residue "+str(i[1])+" and name CA and (resname  GLY PRO)")):
        res2=topology.select("residue "+str(i[1])+" and name CA and (resname  GLY PRO)")
        #print(i[1],res2)
    if (topology.select("residue "+str(i[0])+" and name CB")):
        res1=topology.select('name CB and residue '+str(i[0]))
        #print(i[0],res1)
    if (topology.select("residue "+str(i[1])+" and name CB")):
        res2=topology.select('name CB and residue '+str(i[1]))
        #print(i[1],res2)
    lista_CB.append([int(res1),int(res2)])

print('CB-atomlist')
for j in lista_CB:
    print(j[0],j[1])
print(len(lista_CB))
print(lista_CB)

for i in lista:
    print(i)


      
traj = md.load(xtc, top=pdb, stride=50)
#distance_nm=md.compute_contacts(traj, lista, scheme='ca')
distance_nm=md.compute_distances(traj,lista_CB)

#This is converting to Angstrom
calc=distance_nm[:] * 10
print(calc)

with open(data_f+'AF_contacts_calc.txt', 'w') as f:
    for line in range(len(calc[:,0])):

        f.write(str(line)+" ")
        for j in range(len(lista)): 
            f.write(str(calc[line,j])+" ")
        
        f.write("\n")        
######
#traj[1]
topology = md.load(pdb).topology
print(topology)

table, bonds = topology.to_dataframe()
print(table.head())


# import libraries

i=1
# plot experimental average and error
_ = plt.axvline(exp[i-1,0],c='k',label="AF2")
_ = plt.axvline(exp[i-1,0]-exp[i-1,1],c='k',linestyle="--")
_ = plt.axvline(exp[i-1,0]+exp[i-1,1],c='k',linestyle="--")

# Plot calculated average
average = statistics.median(calc[:,i-1])
_ = plt.axvline(average,c='r',label="Prior ensemble",lw=2)

# plot histogram of the data
_ = plt.hist(calc[:,i-1],bins=100,density=True,alpha=0.3,color='r')
plt.title(labels[i-1])
plt.xlabel("Distance (A)",size=15) 
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.legend()
plt.savefig(data_f+'r_ieq1.pdf')
plt.clf()
plt.close()

i=2    
# plot experimental average and error
_ = plt.axvline(exp[i-1,0],c='k',label="AF2")
_ = plt.axvline(exp[i-1,0]-exp[i-1,1],c='k',linestyle="--")
_ = plt.axvline(exp[i-1,0]+exp[i-1,1],c='k',linestyle="--")

# Plot calculated average
average = statistics.median(calc[:,i-1])
_ = plt.axvline(average,c='r',label="Prior ensemble",lw=2)

# plot histogram of the data
_ = plt.hist(calc[:,i-1],bins=100,density=True,alpha=0.3,color='r')

plt.title(labels[i-1])
plt.legend()

plt.xlabel("Distance (A)",size=15) 
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(data_f+'r_ieq2.pdf')
plt.clf()
plt.close()
##############

bme_dir = os.getcwd().split("notebook")[0]
sys.path.append(bme_dir)
print(bme_dir)
import BME as BME
exp_file = data_f+"AF_contacts_constr.txt"
calc_file = data_f+"AF_contacts_calc.txt" 

print(exp_file)
print(calc_file)

########

# initialize. A name must be specified 

rew = BME.Reweight(data_f+"example_AF_IDP")

print(rew)
# load the experimental and calculated datasets
rew.load(exp_file,calc_file)

#theta_fin=0
thetas = [0.01,0.05,0.1,0.2,0.5,0.7,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29.30,100,1000,10000,100000]

fig, ax = plt.subplots(figsize=(20, 10))

chi2 = []
phis = []
for t in thetas:
    chi2_before, chi2_after, phi = rew.fit(theta=t)
    print(phi)
    phis.append(phi)
    chi2.append(chi2_after/chi2_before)
    #(t,chi2_after/chi2_before,c="k")
    #plt.scatter(t,phi,c="r")
plt.plot(thetas,phis,"-o",label="Phi",c="k")
plt.plot(thetas,chi2,"-o",label="Chi2 reduction",c="r")
plt.legend()
plt.xscale('log')
plt.xlabel("Theta",size=15)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.legend(fontsize=15)
plt.savefig(data_f+'optim_theta.pdf',bbox_inches='tight')
#plt.clf()
plt.close()


chi2 = []
phis = []
for t in thetas:
    chi2_before, chi2_after, phi = rew.fit(theta=t)
    print(phi)
    if (phi>0.05 and phi<0.95):
        print('in if')
        theta_fin=t
        print('THETA fin',theta_fin)
        break

chi2_before, chi2_after, phi = rew.fit(theta=theta_fin)

stats = rew.predict(exp_file,calc_file,data_f+"example_AF_IDP")
f=open(data_f+'opt.log','w')
f.write("%10s %10s %10s" % (" ","Original","Optimized"))
f.write('\n')
f.write("%10s %10.3f %10.3f" % ("Chi2",stats[0],stats[3]))
f.write('\n')
f.write("%10s %10.3f %10.3f" % ("RMSD",stats[1],stats[4]))
f.write('\n')
f.write("%10s %10d %10d" % ("Violations",stats[2],stats[5]))
f.write('\n')
f.write("Fraction of effective frames %.2f" % phi)
f.write('\n')
f.write("Theta fin %.4f" % theta_fin)
f.write('\n')
f.close()
w0 = rew.get_w0()
w_new = rew.get_weights()
print("--------------------------------------")
print("")

########

results = np.loadtxt(data_f+"example_AF_IDP",usecols=(1,2,3,4,5))
labels = np.loadtxt(data_f+"example_AF_IDP",usecols=(0),dtype=str)

# order by magnitude to make the plot nicer
idx_ordered = np.argsort(results[:,0])


fig, ax = plt.subplots(figsize=(20, 10))
xx = range(len(labels))
# plot experiment
plt.errorbar(xx, results[idx_ordered,0],results[idx_ordered,1],c='k',fmt="o",label="_no_legend")
plt.scatter(xx, results[idx_ordered,0],c='k',label="AF")


plt.scatter(xx,results[idx_ordered,2],c='r',label="CG")
plt.scatter(xx,results[idx_ordered,3],c='b',label="AF-CG")
plt.yticks(fontsize=15)
plt.legend(fontsize=15,loc="upper left") 

plt.ylabel("Distance (A)",size=15)
_ = plt.xticks(xx,[labels[l] for l in idx_ordered],rotation=90,size=11)


plt.savefig(data_f+'Prior_Posterior_AF2_distances.pdf',bbox_inches='tight')
plt.clf()
plt.close()
#######

i=2
# plot experimental average and error
_ = plt.axvline(results[i-1,0],c='k',label="AF")
_ = plt.axvline(results[i-1,0]-results[i-1,1],c='k',linestyle="--")
_ = plt.axvline(results[i-1,0]+results[i-1,1],c='k',linestyle="--")

# Plot calculated average
average = statistics.median(calc[:,i-1])
_ = plt.axvline(average,c='r',label="CG",lw=2)

print(len(w_new))
average_optimized = weighted_median(calc[:,i-1],weights=w_new)
_ = plt.axvline(average_optimized,c='b',label="AF-CG",lw=2)

# plot histogram of the data
_ = plt.hist(calc[:,i-1],bins=100,density=True,alpha=0.3,color='r')

_ = plt.hist(calc[:,i-1],bins=100,density=True,alpha=0.3,color='b',weights=w_new)
plt.xlim(0,40)
plt.title(labels[i-1])
plt.xlabel("Distance (A)",size=15) 
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.ylabel('P',fontsize=15)
plt.legend()
plt.savefig(data_f+'r_ieq2_post.pdf',bbox_inches='tight')
plt.clf()
plt.close()


#######
i=2
# plot experimental average and error
_ = plt.axvline(results[i-1,0],c='k',label="AF")
_ = plt.axvline(results[i-1,0]-results[i-1,1],c='k',linestyle="--")
_ = plt.axvline(results[i-1,0]+results[i-1,1],c='k',linestyle="--")

# Plot calculated average
average = statistics.median(calc[:,i-1])
_ = plt.axvline(average,c='r',label="CG",lw=2)

print(len(w_new))
average_optimized = weighted_median(calc[:,i-1],weights=w_new)
_ = plt.axvline(average_optimized,c='b',label="AF-CG",lw=2)

# plot histogram of the data
_ = plt.hist(calc[:,i-1],bins=100,density=True,alpha=0.3,color='r')

_ = plt.hist(calc[:,i-1],bins=100,density=True,alpha=0.3,color='b',weights=w_new)
plt.xlim(0,40)
plt.title(labels[i-1])
plt.xlabel("Distance (A)",size=15) 
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.ylabel('P',fontsize=15)
plt.legend()
plt.savefig(data_f+'r_ieq2_post.pdf',bbox_inches='tight')
plt.clf()
plt.close()
##########RG
#Rg
# Plot calculated average
i=0
average = np.average(rg_A[:])
_ = plt.axvline(average,c='r',label="CG",lw=2)

print(len(w_new))
average_optimized =np.average(rg_A[:],weights=w_new)
_ = plt.axvline(average_optimized,c='b',label="AF-CG",lw=2)

# plot histogram of the data
_ = plt.hist(rg_A[:],bins=100,density=True,alpha=0.3,color='r')

_ = plt.hist(rg_A[:],bins=100,density=True,alpha=0.3,color='b',weights=w_new)


plt.xlabel('Radius of gyration (A)',fontsize=15)
plt.ylabel('P',fontsize=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.legend()


plt.savefig(data_f+'R_g.pdf',bbox_inches='tight')
plt.clf()
plt.close()

#print(average_optimized)
counts, bins = np.histogram(rg_A[:])
mids = 0.5*(bins[1:] + bins[:-1])
probs = counts / np.sum(counts)

mean = np.sum(probs * mids)  
sd = np.sqrt(np.sum(probs * (mids - mean)**2))

print(average,mean,sd)

counts, bins = np.histogram(rg_A[:],weights=w_new)
mids = 0.5*(bins[1:] + bins[:-1])
probs = counts / np.sum(counts)

mean = np.sum(probs * mids)  
sd = np.sqrt(np.sum(probs * (mids - mean)**2))

print(average_optimized,mean,sd)

##########




dssp=md.compute_dssp(traj,simplified='false')
print('dssp data shape', dssp.shape)
residuesALL=b = [i+1 for i in range(nresid)]
array2 = np.vstack((residuesALL,dssp))
dssp=np.c_[["" for x in range(len(array2))],array2]

df2=pd.DataFrame(data=dssp[1:,1:],
                index=dssp[1:,0],
                columns=dssp[0,1:])

print(pd.DataFrame(data=dssp[1:,1:],
                  index=dssp[1:,0],
                  columns=dssp[0,1:]))

df2.to_csv(data_f+'dssp.dat', header=True, index=None, sep=' ', mode='w')


print(len(dssp)-1)

dsspC_resid=[]
dsspH_resid=[]
dsspE_resid=[]

for resid in range(1,nresid+1):
    dummyC=0
    dummyH=0
    dummyE=0
    for t in range(0,len(dssp)-1):
        t_i=t+1
       #print(weights_10[t],sasa[t_i][resid])
        if (dssp[t_i][resid] =='C'):
            dummyC+=float(w0[t])
        elif(dssp[t_i][resid] =='H'):
            dummyH+=float(w0[t])
        elif(dssp[t_i][resid] =='E'):
            dummyE+=float(w0[t])
        
    dsspC_resid.append(dummyC)
    dsspH_resid.append(dummyH)
    dsspE_resid.append(dummyE)
    print(resid,dummyC,dummyH,dummyE)


########
lista_all=[]
table, bonds = topology.to_dataframe()
print(len(table))

for i in range(0,len(table)-1):
    for j in range(i+1,len(table)):
        lista_all.append([i,j])

distance_nm=md.compute_distances(traj,lista_all)
calc_all=distance_nm[:] * 10
#ev 100 works better
plt.plot(ppdf_SAXS_bin[:-1], ppdf_SAXS_prob[:-1]/sum(ppdf_SAXS_prob[:-1]),label="SAXS",c='black')
#plt.scatter(ppdf_SAXS_bin[:-1], ppdf_SAXS_prob[:-1]/sum(ppdf_SAXS_prob[:-1]), marker='+',lw=1,c='black')
yerror=SAXS_err[:-1]
plt.errorbar(ppdf_SAXS_bin[:-1], ppdf_SAXS_prob[:-1]/sum(ppdf_SAXS_prob[:-1]),yerr=yerror/sum(ppdf_SAXS_prob[:-1]),c='black')


flat_all=calc_all.flatten()
histCG_all, bin_edges = np.histogram(flat_all,bins=ppdf_SAXS_bin)

print(len(bin_edges[:-1]),len(histCG_all))
plt.plot(bin_edges[:-1],histCG_all/(sum(histCG_all)),label="CG",c='red')
#plt.scatter(bin_edges[:-1],histCG_all/(sum(histCG_all)), marker='+',lw=1,c='red')


w_new_big=np.array([])
w_new_ev=w_new

for i in range(0,len(w_new_ev)):
    #lst=[w_new_ev10[i]]*calc_all.shape[1]
    w_new_big=np.append(w_new_big,[w_new_ev[i]]*calc_all.shape[1])
print(len(w_new_big))
w_new_big.flatten()
histCG_all_AF, bin_edges = np.histogram(flat_all,bins=ppdf_SAXS_bin,weights=w_new_big)
plt.plot(bin_edges[:-1],histCG_all_AF/sum(histCG_all_AF),label="AF-CG",c='blue')
#plt.scatter(bin_edges[:-1],histCG_all_AF/sum(histCG_all_AF), marker='+',lw=1,c='blue')

plt.xlabel("r ($\AA$)",fontsize=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.ylabel('P',fontsize=15)
plt.legend()
plt.savefig(data_f+'PDDF_CG_AF.pdf',bbox_inches='tight')
plt.clf()
plt.close()


from scipy.special import rel_entr


def chi2_distance(A, B):
 
    # compute the chi-squared distance using above formula
    chi = 0.5 * np.sum([((a - b) ** 2) / (a + b)
                      for (a, b) in zip(A, B)])
 
    return chi
 


def smooth(arr,window_size):
    i = 0
    # Initialize an empty list to store moving averages
    moving_averages = []
    while i < len(arr) - window_size + 1:
    
    # Store elements from i to i+window_size
    # in list to get the current window
        window = arr[i : i + window_size]
  
    # Calculate the average of current window
        window_average = round(sum(window) / window_size, 5)
      
    # Store the average of current
    # window in moving average list
        moving_averages.append(window_average)
      
    # Shift window to right by one position
        i += 1
    return(moving_averages)

def return_intersection(hist_1, hist_2):
    minima = np.minimum(hist_1, hist_2)
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
    return intersection


w_s=2
#Raw
exp = abs(ppdf_SAXS_prob[:-1]/sum(ppdf_SAXS_prob[:-1]))
sim = histCG_all/sum(histCG_all)
exp=[i if i != 0 else 0.0000001 for i in exp]
sim=[i if i != 0 else 0.0000001 for i in sim]
#Smooth
ma_exp=smooth(exp,w_s)
ma_sim=smooth(sim,w_s)

Overlap_CG=return_intersection(sim,exp)
DKL_CG=sum(rel_entr(sim, exp))
#print(sim,exp)


#print(sum(sim),sum(exp))
f=open(data_f+'chi2_stats.dat','w')
####
chi2_CG_SAXS = chi2_distance(exp,sim)
chi2_CG_SAXS_sm = chi2_distance(ma_sim,ma_exp)
f.write("The CG Chi-square distance is : "+ str(chi2_CG_SAXS))
f.write("\n")

plt.plot(ppdf_SAXS_bin[:-1], ppdf_SAXS_prob[:-1]/sum(ppdf_SAXS_prob[:-1]),label="SAXS",c='black')
plt.plot(bin_edges[:-w_s],ma_sim,label="CG smooth",c='red')
#########

#AF 
sim_AF = histCG_all_AF/sum(histCG_all_AF)
sim_AF=[i if i != 0 else 0.0000001 for i in sim_AF]

Overlap_CG_AF=return_intersection(sim_AF,exp)
DKL_CG_AF=sum(rel_entr(sim_AF, exp))

plt.plot(ppdf_SAXS_bin[:-w_s],ma_exp,label="SAXS smooth",c='yellow')

ma_exp=smooth(exp,w_s)
ma_sim_AF=smooth(sim_AF,w_s)






Overlap_CG_AF=return_intersection(sim_AF,exp)
chi2_AF_CG_SAXS = chi2_distance(sim_AF,exp)
f.write("The AF_CG Chi-square distance is: "+str(chi2_AF_CG_SAXS))
f.write("\n")
chi2_AF_CG_SAXS_sm = chi2_distance(ma_sim_AF,ma_exp)
f.write("The CG Chi-square smooth distance is: "+str(chi2_CG_SAXS_sm))
f.write("\n")
f.write("The AF CG Chi-square smooth distance is : "+str(chi2_AF_CG_SAXS_sm))
f.write("\n")
#scipy.stats.chisquare(sim, f_exp=exp)
f.write("The CG Overlap "+str(Overlap_CG))
f.write("\n")
f.write("The AF_Overlap "+str(Overlap_CG_AF))
f.write("\n")

f.write("DKL CG "+str(DKL_CG))
f.write("\n")
f.write("DKL CG-AF "+str(DKL_CG_AF))
f.close()
plt.plot(bin_edges[:-w_s],ma_sim_AF,label="AF-CG smooth",c='blue')
plt.legend()
plt.savefig(data_f+'PDDF_CG_AF_smooth.pdf',bbox_inches='tight')
plt.clf()
plt.close()
######

print(len(dssp)-1)

dsspC_resid_w=[]
dsspH_resid_w=[]
dsspE_resid_w=[]

for resid in range(1,nresid+1):
    dummyC=0
    dummyH=0
    dummyE=0
    for t in range(0,len(dssp)-1):
        t_i=t+1
       #print(weights_10[t],sasa[t_i][resid])
        if (dssp[t_i][resid] =='C'):
            dummyC+=float(w_new[t])
        elif(dssp[t_i][resid] =='H'):
            dummyH+=float(w_new[t])
        elif(dssp[t_i][resid] =='E'):
            dummyE+=float(w_new[t])
        
    dsspC_resid_w.append(dummyC)
    dsspH_resid_w.append(dummyH)
    dsspE_resid_w.append(dummyE)
    print(resid,dummyC,dummyH,dummyE)


###
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(residuesALL,dsspC_resid,color="r",alpha=0.3,label="CG")
ax.bar(residuesALL,dsspC_resid_w,color="b",alpha=0.3,label="AF-CG`")
plt.ylim(0,1.2)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

plt.xlabel('Residue',fontsize=15)
plt.ylabel('% coil',fontsize=15)
#plt.legend(fontsize=20,loc=1,prop={'size': 5})
ax.tick_params(axis='both', labelsize=15)
ax.legend()
plt.savefig(data_f+'coil.pdf',bbox_inches='tight')
plt.clf()
plt.close()



#######
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(residuesALL,dsspH_resid,color="r",alpha=0.3,label="CG")
ax.bar(residuesALL,dsspH_resid_w,color="b",alpha=0.3,label="AF-CG")
plt.ylim(0,1.2)


plt.xlabel('Residue',fontsize=15)
plt.ylabel('% helix',fontsize=15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

#plt.legend(fontsize=20,loc=1,prop={'size': 5})
ax.tick_params(axis='both', labelsize=15)
ax.legend()
plt.savefig(data_f+'helix.pdf',bbox_inches='tight')
plt.clf()
plt.close()


#######
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(residuesALL,dsspE_resid,color="r",alpha=0.3,label="CG")
ax.bar(residuesALL,dsspE_resid_w,color="b",alpha=0.3,label="AF-CG")
plt.ylim(0,1.2)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.xlabel('Residue',fontsize=15)
plt.ylabel('% β-sheet',fontsize=15)
#plt.legend(fontsize=20,loc=1,prop={'size': 5})
ax.tick_params(axis='both', labelsize=15)
ax.legend()
plt.savefig(data_f+'beta.pdf',bbox_inches='tight')
plt.clf()
plt.close()
####
