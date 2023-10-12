import itertools
import os
from pathlib import Path

import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
#from matplotlib.colors import LinearSegmentedColormap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import mdtraj as md
#import seaborn as sns
#from scipy.stats import gaussian_kde

#rc("font", **{"family": "sans-serif",
#colors = sns.husl_palette()

#get_ipython().run_line_magic('matplotlib', 'inline')
#get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")


# In[5]:


#from typing import Union, List, Callable


KBT = 2.494339 # kJ/mol, 300K
n_replica = 1
top_file = "input_af.pdb"
biased_cvs = {}
time=np.loadtxt("FULLBIAS")[:, 0]
bias=np.loadtxt("FULLBIAS")[:, 1]
bias10=bias
print(len(bias10))
weights = np.exp(bias10/ KBT)
weights /= weights.sum()
nframes = {weights.shape[0]}


# ## Clustering Full Protein
# Gromos clustering with a trajectory subsample, based on CA RMSD. We first evaluate different cutoffs:
# gmx trjconv -f cat_traj_recon.gro -s structure.pdb -o cat_traj_recon_noh_skip10_noh.gro -skip 10
#  

# In[15]:


#traj = md.load_xtc("conf_0_recon.xtc",top='output.pdb')
traj = md.load("conf_0_recon.gro",stride=1)


# In[25]:


N = 1000
n_rounds = 5

folder = Path(f"./")
n_frames = traj.n_frames
inds = np.arange(n_frames)

for i in range(n_rounds):
    (folder / f"r{i}").mkdir(exist_ok=True)
    inds_sample = np.random.choice(inds, size=N, replace=False,
                               p=weights / weights.sum())
    traj[inds_sample].save_xtc((folder / f"r{i}" / f"conf-protein.xtc").as_posix())
