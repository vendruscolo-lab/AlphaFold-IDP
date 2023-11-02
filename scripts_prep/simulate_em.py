from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import pandas as pd
import numpy as np
import pickle
import mdtraj as md
from mpi4py import MPI
from openmmplumed import PlumedForce
comm1 = MPI.COMM_SELF
comm2 = MPI.COMM_WORLD

platform = Platform.getPlatformByName('CPU')
#properties={}

residues = pd.read_csv("residues.csv").set_index('one')

fasta_file = open('sequence.dat',mode='r')
fasta = fasta_file.read().replace("\n","")
fasta_file.close()

#fasta = """
#MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGAGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPAKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRACDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSAGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGAGSASNAGSGSGFNGGFGSSMDSKSSGAGM""".replace('\n', '')

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

residues = adjust_terminals_HIS(residues,7.4,fasta)

atomic_number = 117
for i,r in residues.iterrows():
    name = 'carbon-{}'.format(i)
    symbol = 'C-{}'.format(i)
    new_element = Element(atomic_number,name,symbol,r.MW*dalton)
    atomic_number += 1

pdb = PDBFile('input_af.pdb')

top = pdb.topology

for chain in top._chains:
    for residue in chain._residues:
        for atom in residue._atoms:
            atom.element = Element._elements_by_symbol['C-{}'.format(residue.name)]
            
atoms = list(top.atoms())

for i in range(len(fasta)-1):
    top.addBond(atoms[i],atoms[i+1])
    
forcefield = ForceField("forcefield.xml")

system = forcefield.createSystem(top,nonbondedMethod=CutoffNonPeriodic)

system.removeForce(3)


with open('r1_excl.pkl', 'rb') as fp:
    r1_exclusions = pickle.load(fp)
    
forces = system.getForces()

forces[1].setCutoffDistance(4*nanometers)
for bond in r1_exclusions:
    forces[1].addExclusion(bond[0],bond[1])
    
forces[2].setForceGroup(1)
forces[2].setCutoffDistance(2*nanometers)
for bond in r1_exclusions:
    forces[2].addExclusion(bond[0],bond[1])

#open text file in read mode
                        
N_res = len(fasta)
N_save = 3000 if N_res < 100 else int(np.ceil(3e-4*N_res**2)*1000)
N_steps = 1100*N_save
DCD_output_file = "output.dcd"
stats_output_file = "stats.csv"

integrator = LangevinIntegrator(298*kelvin, 0.01/picosecond, 0.005*picoseconds)


#simulation = Simulation(top, system, integrator,platform,properties)
simulation = Simulation(top, system, integrator,platform)
    
simulation.context.setPositions(pdb.positions)

simulation.minimizeEnergy(0)

simulation.saveCheckpoint("checkpoint")
