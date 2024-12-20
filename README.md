# AlphaFold-IDP
This repository provides with custom code and analysis scripts to generate structural ensembles of intrinsically disordered proteins by using AlphaFold generated residue distance maps as constraints using prior-structural ensembles. AlphaFold structural ensemble of Intrinsically Disordered Proteins (IDPs) and Partially Disordered Proteins (PDPs) and figures generated in the publication can be found in [Zenodo](https://zenodo.org/records/14387133). The notebook in `prep_run/prepare_run.ipynb` is an example that generates the AlphaFold MetaInference (AF-MI) structural ensemble of TDP-43 as well provides with free energy surface analysis. For this example, the AF_DATA folder contains AlphaFold predictions generated by using as input the TDP-43WtoA sequence. The expected output can be found in [Zenodo](https://zenodo.org/records/14387133) repository.

## Layout
- `AlphaFold prediction of pairwise distances`: In [this](https://github.com/zshengyu14/ColabFold_distmats/blob/main/AlphaFold2.ipynb) google colab notebook, the user can generate AlphaFold predictions  in terms of means, standard deviation and even of probability distribution of individual pairwise aminoacid distances, from arbitrary aminoacid sequences. The output folder of this procedure is stored in the `AF_DATA` folder.
- `prep_run`: This folder contains the necessary scripts to a) convert atomistic AF pdb to coarse grained pdb, b) create the necessary plumed file in order to run AF-MI, c) run AF-MI by using our new [OPENMM-PLUMED-MPI implementation](https://github.com/vendruscolo-lab/OpenMM-Plumed-MPI) of CALVADOS coarse grained model, d) reconstruct the coarse grained structural ensemble of the disordered protein, f) backmap it to atomistic models e) calculate free energy profiles along predefined collective variables


## Reproducability information 

The analysis was performed on a single machine containing 12 Intel(R) Xeon(R) W-2133 CPU @ 3.60GH CPUs, with total system memory of 100 GB 
Running the notebooks that performs the AF-ensemble generation takes approximately < 10'
The google-colab notebook a few minutes to generate the the AF pairwise aminoacid distance prediction of ~ 100 aminoacid  long proteins.
The backmap analysis took up to to 3 hrs depending on the proteins tested.

## Software requirements

The code can be used in Linux or macOS. 

### Dependencies: 
```
numpy
pandas
scipy
notebook
mdtraj 
matplolib
plumed
gromacs
```

### Installation:

The user should download this github repo locally by `gitclone https://github.com/vendruscolo-lab/AlphaFold-IDP.git`

Then, install and activate an conda environment with [OPENMM-PLUMED-MPI implementation](https://github.com/vendruscolo-lab/OpenMM-Plumed-MPI)

Then, update that environment with the dependencies listed above

## Authors
[Faidon Brotzakis (@fbrotzakis)](https://github.com/fbrotzakis)

[Shengyu Zhang (@zshengyu14)](https://github.com/zshengyu14)

[Mhd. Hussein Murtada (@husseinmur)](https://github.com/husseinmur)

[Michele Vendruscolo (@vendruscolo-lab)](https://github.com/vendruscolo-lab)
## Article

https://www.biorxiv.org/content/10.1101/2024.11.09.622758v1
