#python gen_xml_and_constraints.py sequence.dat 7.4 298 0.2
import pandas as pd
import numpy as np
import pickle
import mdtraj as md
import sys

def adjust_terminals_HIS(r,pH,fasta):
    r.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['X','q'] = r.loc[fasta[0],'q'] + 1.
    r.loc['X','MW'] = r.loc[fasta[0],'MW'] + 2.
    r.loc['X','three'] = 'X'
    
    r.loc['Z'] = r.loc[fasta[-1]]
    r.loc['Z','q'] = r.loc[fasta[-1],'q'] - 1.
    r.loc['Z','MW'] = r.loc[fasta[-1],'MW'] + 16.
    r.loc['Z','three'] = 'Z'
    
    return r.set_index('three')
	
def calculate_yukawa_params(temp,ionic):
    RT = 8.3145*temp*1e-3
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/RT
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    return yukawa_kappa, RT, lB
	
def generate_atomtypes_block(res):
    tag = 'AtomTypes'
    atomtypes = ""
    for i,r in res.iterrows():
        atomtype = '<Type name="{}" class="C" element="C-{}" mass="{}"/>\n'.format(i,i,str(r.MW))
        atomtypes += "  " + atomtype
    out = " <{0}>\n{1} </{0}>\n\n".format(tag,atomtypes)
    return out
	
def generate_residues_block(res):
    tag = 'Residues'
    residues = ""
    for i,r in res.iterrows():
        atom = '   <Atom name="{0}" type="{0}" charge="{1}"/>'.format(i,r.q)
        bond = '   <ExternalBond atomName="{}"/>'.format(i)
        if i not in ['X','Z']:
            bond += '\n   <ExternalBond atomName="{}"/>'.format(i)
        residue = '<Residue name="{}">\n{}\n{}\n  </Residue>\n'.format(i,atom,bond)
        residues += "  " + residue
    out = " <{0}>\n{1} </{0}>\n\n".format(tag,residues)
    return out
	
def generate_harmonicBond_block(length,k):
    tag = "HarmonicBondForce"
    bond = '  <Bond class1="C" class2="C" length="{}" k="{}"/>'.format(length,k)
    out = " <{0}>\n{1}\n </{0}>\n\n".format(tag,bond)
    return out
	
def generate_yukawa_block(res,temp,ionic,r_cut):
    kappa, RT, lB = calculate_yukawa_params(temp,ionic)
    energy_eq = "RT*lB*epsilon1*epsilon2*(exp(-r*kappa)/r-exp(-r_cut*kappa)/r_cut)"
    tag = "CustomNonbondedForce"
    outer = ' <{} energy="{}" bondCutoff="0">\n'.format(tag,energy_eq)
    inner = '  <GlobalParameter name="RT" defaultValue="{}"/>\n'.format(RT)
    inner += '  <GlobalParameter name="lB" defaultValue="{}"/>\n'.format(lB)
    inner += '  <GlobalParameter name="kappa" defaultValue="{}"/>\n'.format(kappa)
    inner += '  <GlobalParameter name="r_cut" defaultValue="{}"/>\n'.format(r_cut)
    inner += '  <PerParticleParameter name="epsilon"/>\n'
    for i,r in res.iterrows():
        inner += '  <Atom type="{}" epsilon="{}"/>\n'.format(i,r.q)
    out = outer+inner+" </{}>\n\n".format(tag)
    return out

def generate_ashbaugh_block(res,rcut):
    lj_eps = 4.184*.2
    threshold= 2**(1/6)
    energy_eq = "c1*vlj+c2-shift;vlj=4*epsilon*((sigma/r)^12-(sigma/r)^6);c1=select(delta(d),lambda,s1);c2=select(delta(d),0,s2);s1=select(step(d),lambda,1);s2=select(step(d),0,(1-lambda)*epsilon);d=r-threshold*sigma;shift=lambda*vlj_shift;vlj_shift=4*epsilon*((sigma/rcut)^12-(sigma/rcut)^6);lambda=0.5*(lambda1+lambda2);sigma=0.5*(sigma1+sigma2)"
    tag = "CustomNonbondedForce"
    outer = '  <{} energy="{}" bondCutoff="0">\n'.format(tag,energy_eq)
    inner = '  <GlobalParameter name="epsilon" defaultValue="{}"/>\n'.format(lj_eps)
    inner += '  <GlobalParameter name="threshold" defaultValue="{}"/>\n'.format(threshold)
    inner += '  <GlobalParameter name="rcut" defaultValue="{}"/>\n'.format(rcut)
    inner += '  <PerParticleParameter name="sigma"/>\n'
    inner += '  <PerParticleParameter name="lambda"/>\n'
    for i,r in res.iterrows():
        inner += '  <Atom type="{}" sigma="{}" lambda="{}"/>\n'.format(i,r.sigmas,r.lambdas)
    out = outer+inner+" </{}>\n\n".format(tag)
    return out

def generate_forcefield_xml(fasta, residues, xml_name, pH, temp, ionic):
    residues = adjust_terminals_HIS(residues,pH,fasta)
    atomTypes_block = generate_atomtypes_block(residues)
    residues_block = generate_residues_block(residues)
    harmonic_bond_block = generate_harmonicBond_block(0.38,8033.0)
    yukawa_block = generate_yukawa_block(residues,temp,ionic,4)
    ashbaugh_block = generate_ashbaugh_block(residues,2)
    xml_output="<ForceField>\n\n{}{}{}{}{}</ForceField>".format(atomTypes_block,residues_block,harmonic_bond_block,yukawa_block,ashbaugh_block)
    with open(xml_name, "w") as f:
        f.write(xml_output)
		
def create_exclusions_file(n_res):
    excl = []
    for i in range(n_res-1):
        excl.append([i,i+1])
    with open('r1_excl.pkl', 'wb') as fp:
        pickle.dump(excl, fp)

if len(sys.argv) == 5:
    pH = float(sys.argv[2])
    #print(pH)
    temp = float(sys.argv[3])
    ionic = float(sys.argv[4])
else:
    pH = 7.4
    temp = 298
    ionic = 0.2

fasta_file = open(sys.argv[1],mode='r')
fasta = fasta_file.read().replace("\n","")
fasta_file.close()

residues = pd.read_csv('residues.csv').set_index('one',drop=False)

generate_forcefield_xml(fasta, residues, "forcefield.xml", pH, temp, ionic)
create_exclusions_file(len(fasta))
